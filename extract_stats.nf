// Folder containing the input files
params.input_fld = "results/test_dataset"

// input files directory
input_dir = file(params.input_fld)
assert input_dir.isDirectory()

// source scripts directory
script_dir = file("$baseDir/scripts")


// reads input channel. Has items [n. vial, timepoint, file]
regex = /\/vial_(\d+)/
input_ch = Channel.fromPath("${input_dir}/vial_*", type:"dir")
    .map { it -> [vial: (it =~ regex)[0][1], fld: it] }

scripts = ["extract_consensus_freq.py", "extract_gap_freq.py"]


process extract_script {
    label 'q6h_1core'

    publishDir "results/${input_dir.getName()}/vial_${vial}/stats/", mode: 'copy'

    input:
        tuple val(vial), path(vial_fld) from input_ch
        each script from scripts

    output:
        path("stats_table_*.pkl.gz")

    script:
        """
        python3 $script_dir/$script --vial_fld $vial_fld --output_fld . --verbose
        """
}
