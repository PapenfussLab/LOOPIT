process GET_SEGMENTS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(breakpoints), path(reference_regions)

    output:
    tuple val(meta), path("${meta.id}_segments.bed"), emit: segments
    path "versions.yml",                               emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_segments.R \\
        --sample             "${meta.id}" \\
        --breakpoints        "${breakpoints}" \\
        --reference-regions  "${reference_regions}" \\
        --chromosomes        "${params.chromosomes}" \\
        --outfile            "${prefix}_segments.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
