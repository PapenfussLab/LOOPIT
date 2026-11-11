process GET_BREAKPOINTS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(repeat_bed)

    output:
    tuple val(meta), path("${meta.id}_breakpoints.rds"), emit: breakpoints
    path "versions.yml",                                  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode   = meta.mode ?: 'SR'
    """
    get_breakpoints.R \\
        --sample      "${meta.id}" \\
        --vcf         "${vcf}" \\
        --mode        "${mode}" \\
        --repeat-bed  "${repeat_bed}" \\
        --outdir      "." \\
        --outfile     "${prefix}_breakpoints.rds"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
