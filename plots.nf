// Folder containing the input files
params.input_fld = "results/test_dataset"

// input files directory
input_dir = file(params.input_fld.replaceAll('/$',''))
assert input_dir.isDirectory()

// function to extract vial number from folder name
def extract_vial_n(fld) {
    regex = /\/vial_(\d+)/
    match = (fld =~ regex)[0]
    return [vial: match[1] as Integer, fld: fld]
}

// process to perform plots
process plot_script {

    label 'q30m_1core'
    conda 'conda_envs/bioinfo_raw.yml'

    publishDir "figures/${input_dir.getName()}/vial_${vial}/", mode: 'copy'

    input:
        tuple val(vial), path(vial_fld)
        each script

    output:
        path("*.pdf")
        path("*.csv") optional true

    script:
        """
        python3 $baseDir/scripts/$script --vial_fld $vial_fld --fig_fld .
        """
}

workflow {

    // channel for vial folders [vial_n, vial_fld]
    vials = Channel.fromPath("${input_dir}/vial_*", type:"dir")
        .map { it -> extract_vial_n(it) }

    // list of plot scripts
    scripts = [
        "plot_coverage.py",
        "plot_consensus_freq.py",
        "plot_gaps.py",
        "plot_insertions.py",
        "plot_clips.py"
        ]

    // execute plot scripts on each vial
    plot_script(vials, scripts)
}
