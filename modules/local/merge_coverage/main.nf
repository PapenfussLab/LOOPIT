process MERGE_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(coverage_files)

    output:
    tuple val(meta), path("${meta.id}_coverage.rds"), emit: coverage
    path "${meta.id}_merge_coverage_debug.txt",        emit: debug
    path "versions.yml",                               emit: versions

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def cov_list   = coverage_files instanceof List
        ? coverage_files.join(' ')
        : coverage_files.toString()
    """
    merge_coverage.R \\
        --sample  "${meta.id}" \\
        --infiles "${cov_list}" \\
        --outfile "${prefix}_coverage.rds"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}