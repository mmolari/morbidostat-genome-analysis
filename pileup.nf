// Folder containing the input files
params.input_fld = "test_dataset"

// label that marks the first timepoint in the series, from which the reference genome should be taken
// TODO: convert time and vial numbers to integer? This will make it easier to find the first timepoint? 
params.time_beg = "1"

// Parameters for the pileup script
params.qual_min = 15
params.max_insertion_size = 1000

// input files directory
input_dir = file(params.input_fld)
assert input_dir.isDirectory()

// source scripts directory
script_dir = file("$baseDir/scripts")

// string corresponding to the first timepoint
time_beg_id = params.time_beg as String

// Function to extract vial number and timepoint from the path of the input file
def extract_vial_timepoint(x, ending) {
    m = x =~ /\/vial_(\d+)\/time_([^\/\s]+)\/${ending}$/
    return [vial: m[0][1], timepoint: m[0][2], file: x]
}

// reads input channel. Has items [n. vial, timepoint, file]
reads_in = Channel.fromPath("${input_dir}/vial_*/time_*/reads.fastq.gz")
    .map { extract_vial_timepoint(it, "reads.fastq.gz") }

// assembled genome input channel. Has items [n. vial, timepoint, file]
// filtered so that only the first timepoint is kept
genome_in = Channel.fromPath("${input_dir}/vial_*/time_*/assembled_genome/assembly.fna")
    .map { extract_vial_timepoint(it, "assembled_genome/assembly.fna")}
    .filter { it.timepoint == time_beg_id} 


//  combine each set of reads with the input genome from the first timepoint
minimap_in = genome_in.cross(reads_in)
    .map {
        // returns: [vial n., timepoint, reference genome, reads]
        it -> [it[0].vial, it[1].timepoint, it[0].file, it[1].file]
    }

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

    script:
        """
        python3 $script_dir/create_allele_counts.py \
            --bam_file reads.sorted.bam \
            --out_dir pileup \
            --qual_min $params.qual_min \
            --max_insertionsize $params.max_insertion_size
        """
}

