// Folder containing the input files
params.input_fld = "test_dataset"

// label that marks the first timepoint in the series, from which the reference genome should be taken
// TODO: convert time and vial numbers to integer? This will make it easier to find the first timepoint? 
params.time_beg = 1

// Parameters for the pileup script
params.qual_min = 5
params.clip_minL = 100

// input files directory
input_dir = file(params.input_fld)
assert input_dir.isDirectory()

// results directory
output_dir = "$baseDir/results/${input_dir.getName()}"

// function to extract vial and timepoint from file path
def extract_vial_timepoint(fname) {
    regex = /\/vial_(\d+)\/time_(\d+)\//
    match = (fname =~ regex)[0]
    return [vial: match[1] as Integer, timepoint: match[2] as Integer]
}

// process to create symlink of reference genome in each vial/timepoint folder
process create_symlinks {

    publishDir "$output_dir/vial_${vial}/time_${timepoint}/",
        mode: "copyNoFollow"

    input:
        tuple val(vial), val(timepoint), path("ref_genome.fa"), path("ref_genome.gbk"), path("reads.fastq.gz")

    output:
        path "ref_genome.fa"
        path "ref_genome.gbk"
        path "reads.fastq.gz"

    """
    echo "creating symlink for reference genome for vial $vial timepoint $timepoint"
    """
}

// map reads to the reference genome, the output is a sam file
process map_reads {

    label 'q30m'
    conda "conda_envs/read_map.yml"

    input:
        tuple val(vial), val(timepoint), path(genome_fa), path(genome_gbk), path(reads)

    output:
        tuple val(vial), val(timepoint), path("reads.sam")

    script:
        """
        minimap2 -a -x map-ont -t ${task.cpus} $genome_fa $reads > reads.sam
        """
}


// transform the sam file into a sorted bam file
process sort_mapped_reads {

    label 'q30m'
    conda "conda_envs/read_map.yml"

    input:
        tuple val(vial), val(timepoint), path("reads.sam")

    output:
        tuple val(vial), val(timepoint), path("reads.sorted.bam")

    script:
        """
        samtools sort -@ ${task.cpus} reads.sam > reads.sorted.bam
        """
}

// create an index for the bam file
process index_sorted_reads {

    label 'q30m'
    conda "conda_envs/read_map.yml"


    publishDir "$output_dir/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), path("reads.sorted.bam")

    output:
        tuple val(vial), val(timepoint), path("reads.sorted.bam"), path("reads.sorted.bam.bai")

    script:
        """
        samtools index -@ ${task.cpus} reads.sorted.bam
        """
}

// Perform the pileup with richard's script
process pileup {

    label 'q6h_1core'
    conda "conda_envs/bioinfo_raw.yml"

    publishDir "$output_dir/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), path("reads.sorted.bam"), path("reads.sorted.bam.bai")


    output:
        path("pileup/allele_counts.npz")
        path("pileup/insertions.pkl.gz")
        path("pileup/clips.pkl.gz")

    script:
        """
        python3 $baseDir/scripts/create_allele_counts.py \
            --bam_file reads.sorted.bam \
            --out_dir pileup \
            --qual_min $params.qual_min \
            --clip_minL $params.clip_minL
        """
}

process unmapped {
    label 'q30m'
    conda 'conda_envs/bioinfo_raw.yml'

    publishDir "$output_dir/vial_${vial}/time_${timepoint}/pileup", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), path("reads.sorted.bam"), path("reads.sorted.bam.bai")


    output:
        path("unmapped.csv"), optional: true
        path("unmapped.fastq.gz"), optional: true


    script:
        """
        python3 $baseDir/scripts/pileup_unmapped.py \
            --bam reads.sorted.bam \
            --df_out unmapped.csv \
            --fastq_out unmapped.fastq.gz
        """
}

process non_primary {
    label 'q30m'
    conda 'conda_envs/bioinfo_raw.yml'

    publishDir "$output_dir/vial_${vial}/time_${timepoint}/pileup", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), path("reads.sorted.bam"), path("reads.sorted.bam.bai")


    output:
        tuple val(vial), val(timepoint), path("non_primary.csv"), optional: true

    script:
        """
        python3 $baseDir/scripts/pileup_secondary_supplementary.py \
            --bam reads.sorted.bam \
            --df_out non_primary.csv \
        """
}


process plot_non_primary_single {
    label 'q30m'
    conda 'conda_envs/bioinfo_raw.yml'

    publishDir "figures/${input_dir.getName()}/vial_${vial}/non_primary", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), path("non_primary.csv")

    output:
        path("*.pdf")

    script:
        """
        python3 $baseDir/scripts/plot_secondary.py --df non_primary.csv \
            --pdf secondary_t_${timepoint}.pdf
        python3 $baseDir/scripts/plot_supplementary.py --df non_primary.csv \
            --pdf supplementary_t_${timepoint}.pdf
        """
}

process plot_non_primary_vs_time {
    label 'q30m'
    conda 'conda_envs/bioinfo_raw.yml'

    publishDir "figures/${input_dir.getName()}/vial_${vial}/non_primary", mode: 'copy'

    input:
        tuple val(vial), val(timepoints), path("non_primary_*.csv")

    output:
        path("*.pdf")

    script:
        """
        python3 $baseDir/scripts/plot_non_primary_vs_t.py \
            --dfs non_primary_*.csv --ts ${timepoints.join(" ")} \
            --pdf_sec secondary_vs_t.pdf --pdf_suppl supplementary_vs_t.pdf
        """
}

workflow pileup_workflow {
    main:

        // string corresponding to the first timepoint
        time_beg_id = params.time_beg as Integer

        // reads input channel. Has items [vial, timepoint, reads]
        reads = Channel.fromPath("${input_dir}/vial_*/time_*/reads.fastq.gz")
            .map {it -> extract_vial_timepoint(it) + [reads:it]}
        
        // assembled genome input channel. Has items [vial, timepoint, fa, gbk]
        // filtered so that only the first timepoint is kept
        assembled_genomes = Channel.fromFilePairs("${input_dir}/vial_*/time_*/assembled_genome/*.{fna,gbk}")
            .map {it -> it[1] } // only file pair
            .map {it -> extract_vial_timepoint(it[0]) + [fna:it[0], gbk:it[1]] }
            .filter {it -> it.timepoint == time_beg_id}


        // combine genomes with reads, using vial as common key.
        // Each set of reads is assigned the rescpective reference genome
        combined = assembled_genomes.cross(reads) {it -> it.vial }
            .map {it -> [it[0].vial, it[1].timepoint, it[0].fna, it[0].gbk, it[1].reads]}

        // create symlinks of reference genomes and reads
        create_symlinks(combined)
        
        // map and sort reads
        sorted_reads = map_reads(combined) | sort_mapped_reads

        // create index for sorted reads
        indexed_reads = index_sorted_reads(sorted_reads)

        // perform pileup and list unmapped and non-primary reads
        pileup(indexed_reads)
        unmapped(indexed_reads)
        df_non_primary = non_primary(indexed_reads)

    emit:
        non_primary = df_non_primary
}

workflow plots_workflow {
    take:
        non_primary
    main:
        plot_non_primary_single(non_primary)
        // [2, [2, 1, 5], [non_primary.csv, non_primary.csv, non_primary.csv]]
        plot_non_primary_vs_time(non_primary.groupTuple())
}

workflow {
    pp_out = pileup_workflow()
    plots_workflow(pp_out.non_primary)
}