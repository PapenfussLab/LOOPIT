process GET_REFERENCE_REGIONS {
    tag "${params.genome}"
    label 'process_medium'
    storeDir "${params.outdir}/reference_data/${params.genome}"

    input:
    path reference_dir

    output:
    path "reference_regions.rds", emit: reference_regions
    path "versions.yml",          emit: versions

    script:
    def fasta_name = params.genome == "GRCh37" ? "Homo_sapiens_assembly19.fasta" : "Homo_sapiens_assembly38.fasta"
    """
    get_reference_regions.R \\
        --reference     "${reference_dir}/${fasta_name}" \\
        --reference-dir "${reference_dir}" \\
        --outfile       "reference_regions.rds"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}