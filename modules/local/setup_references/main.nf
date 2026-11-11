process SETUP_REFERENCES {
    tag "${params.genome}"
    label 'process_single'
    storeDir "${params.outdir}/reference_data/${params.genome}"

    output:
    path "${params.genome}", emit: reference_dir
    path "versions.yml",     emit: versions

    script:
    """
    bash ${projectDir}/scripts/download_reference_data.sh ${params.genome} ${params.genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //' | sed 's/ .*//')
    END_VERSIONS
    """
}