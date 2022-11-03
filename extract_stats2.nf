// folder containing the input files
params.input_fld = "results/test_dataset"

// input files directory
input_dir = file(params.input_fld)
assert input_dir.isDirectory()

// function to extract vial number from folder name
def extract_vial_n(fld) {
    regex = /\/vial_(\d+)/
    match = (fld =~ regex)[0]
    return [vial: match[1] as Integer, fld: fld]
}

// run script to extract info
process extract_script {

    label 'q6h_1core'
    conda 'conda_envs/bioinfo_raw.yml'

    publishDir "${input_dir}/vial_${vial}/stats/", mode: 'copy'

    input:
        tuple val(vial), path(vial_fld)
        each script

    output:
        path("stats_table_*.pkl.gz")

    script:
        """
        python3 $baseDir/scripts/$script --vial_fld $vial_fld --output_fld . --verbose
        """
}

workflow {

    // Channel for vial folders [vial_n, vial_fld]
    vials = Channel.fromPath("${input_dir}/vial_*", type:"dir")
        .map { it -> extract_vial_n(it) }
    
    // set of scripts to run on each vial dataset
    scripts = ["extract_consensus_freq.py", "extract_gap_freq.py"]

    // Execute each script on each vial
    extract_script(vials, scripts)

}