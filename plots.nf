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

// plot_scripts = ["plot_coverage.py", "plot_consensus_frequency.py", "plot_gaps.py", "plot_insertions.py"]
plot_scripts = ["plot_coverage.py"]


process plot_script {
    label 'q30m_1core'

    publishDir "figures/${input_dir.getName()}/vial_${vial}/", mode: 'copy'

    input:
        tuple val(vial), path(vial_fld) from input_ch
        each script from plot_scripts

    output:
        path("*.pdf")
        path("*.csv") optional true

    script:
        """
        python3 $script_dir/$script --vial_fld $vial_fld --fig_fld .
        """
}
