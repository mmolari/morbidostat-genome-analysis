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

    publishDir "$output_dir/vial_${vial}/time_${timepoint}/", mode: 'copy'

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
        path("reads.sorted.bam.bai")

    script:
        """
        samtools index -@ ${task.cpus} reads.sorted.bam
        """
}

// Perform the pileup with richard's script
process pileup {

    label 'q6h_1core'
    conda "conda_envs/bioinfo.yml"

    publishDir "$output_dir/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), path("reads.sorted.bam")

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

workflow {

    // string corresponding to the first timepoint
    time_beg_id = params.time_beg as Integer

    // reads input channel. Has items [vial, timepoint, reads]
    reads = Channel.fromPath("${input_dir}/vial_*/time_*/reads.fastq.gz")
        .map {it -> extract_vial_timepoint(it) + [reads:it]}
    
    // assembled genome input channel. Has items [vial, timepoint, fa, gbk]
    // filtered so that only the first timepoint is kept
    assembled_genomes = Channel.fromFilePairs("${input_dir}/vial_*/time_*/assembled_genome/assembly.{fna,gbk}")
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
    index_sorted_reads(sorted_reads)

    // perform pileup
    pileup(sorted_reads)
}