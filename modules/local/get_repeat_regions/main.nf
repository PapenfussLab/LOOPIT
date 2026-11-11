process GET_REPEAT_REGIONS {
    tag "${params.genome}"
    label 'process_low'
    storeDir "${params.outdir}/reference_data/${params.genome}"

    input:
    path reference_dir

    output:
    path "repeat_mask.bed", emit: repeat_bed
    path "versions.yml",    emit: versions

    script:
    """
    get_repeat_regions.R \\
        --reference-dir "${reference_dir}" \\
        --outfile       "repeat_mask.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
