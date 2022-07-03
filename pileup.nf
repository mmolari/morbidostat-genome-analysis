// Folder containing the input files
params.input_fld = "test_dataset"

// label that marks the first timepoint in the series, from which the reference genome should be taken
// TODO: convert time and vial numbers to integer? This will make it easier to find the first timepoint? 
params.time_beg = "1"

// Parameters for the pileup script
params.qual_min = 5
params.clip_minL = 100

// input files directory
input_dir = file(params.input_fld)
assert input_dir.isDirectory()

// source scripts directory
script_dir = file("$baseDir/scripts")

// string corresponding to the first timepoint
time_beg_id = params.time_beg as String

// reads input channel. Has items [n. vial, timepoint, file]
regex1 = /\/vial_(\d+)\/time_([^\/\s]+)\/reads\.fastq\.gz$/
reads_in = Channel.fromPath("${input_dir}/vial_*/time_*/reads.fastq.gz")
    // .map { extract_vial_timepoint(it) }
    .map {it -> [vial: (it =~ regex1)[0][1],
                timepoint: (it =~ regex1)[0][2],
                file: it]}


// assembled genome input channel. Has items [n. vial, timepoint, file]
// filtered so that only the first timepoint is kept
regex2 = /\/vial_(\d+)\/time_([^\/\s]+)\/assembled_genome\/assembly\.fna$/
genome_in = Channel.fromPath("${input_dir}/vial_*/time_*/assembled_genome/assembly.fna")
    .map {it -> [vial: (it =~ regex2)[0][1],
                timepoint: (it =~ regex2)[0][2],
                file: it]}


//  combine each set of reads with the input genome from the first timepoint
genome_in.cross(reads_in) { it -> it.vial }
    .map {
        // returns: [vial n., timepoint, reference genome, reads]
        it -> [it[0].vial, it[1].timepoint, it[0].file, it[1].file]
    }
    .into { minimap_in; symlink_creation}

// function to create symlinks of reference genomes
def save_genome_symlink(vial, timepoint, genome_fa_file) {

    // make destination folder 
    dest_str = "$baseDir/results/${input_dir.getName()}/vial_${vial}/time_${timepoint}/"
    dest_dir = file(dest_str)
    dest_dir.mkdirs()

    // define names of destination files
    gf_fa_dest = dest_str + "ref_genome.fa"
    gf_gbk_dest = dest_str + "ref_genome.gbk"

    // define genbank source file name
    genome_gbk_file = (genome_fa_file as String).replaceFirst(/assembly.fna$/, "assembly.gbk")

    // create symlinks
    genome_fa_file.mklink(gf_fa_dest, overwrite:true)
    file(genome_gbk_file).mklink(gf_gbk_dest, overwrite:true)

    return "created symlink for reference genome for vial $vial timepoint $timepoint"
}

// create symlinks for reference genome
symlink_creation.map { save_genome_symlink(it[0], it[1], it[2]) } .view()

// map reads to the reference genome, the output is a sam file
process map_reads_to_genome {

    label 'q30m'

    input:
        tuple val(vial), val(timepoint), file(genome), file(reads) from minimap_in 

    output:
        tuple val(vial), val(timepoint), file("reads.sam") into minimap_out

    script:
        """
        minimap2 -a -x map-ont -t ${task.cpus} $genome $reads > reads.sam
        """
}


// transform the sam file into a sorted bam file
process sort_mapped_reads {

    label 'q30m'

    publishDir "results/${input_dir.getName()}/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), file("reads.sam") from minimap_out

    output:
        tuple val(vial), val(timepoint), file("reads.sorted.bam") into samtools_out

    script:
        """
        samtools sort -@ ${task.cpus} reads.sam > reads.sorted.bam
        """
}

// duplicate the channel for the indexing and the pileup
samtools_out.into { index_in; pileup_in}

// Create an index for the bam file
process index_sorted_reads {

    label 'q30m'

    publishDir "results/${input_dir.getName()}/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), file("reads.sorted.bam") from index_in

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

    publishDir "results/${input_dir.getName()}/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        tuple val(vial), val(timepoint), file("reads.sorted.bam") from pileup_in

    output:
        path("pileup/allele_counts.npz")
        path("pileup/insertions.pkl.gz")
        path("pileup/clips.pkl.gz")

    script:
        """
        python3 $script_dir/create_allele_counts.py \
            --bam_file reads.sorted.bam \
            --out_dir pileup \
            --qual_min $params.qual_min \
            --clip_minL $params.clip_minL
        """
}

