process COMPUTE_COVERAGE {
    tag "${meta.id}_${meta.chrom}"
    label 'process_high'

    input:
    tuple val(meta), path(tumour_bed, stageAs: "tumour/*"), path(normal_bed, stageAs: "normal/*"),
          path(segment_regions), path(blacklist), path(repeat_bed)

    output:
    tuple val(meta), path("${meta.id}_${meta.chrom}_coverage.rds"), emit: coverage
    path "${meta.id}_${meta.chrom}_has_normal_debug.txt",            emit: debug
    path "versions.yml",                                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.chrom}"
    """
    get_coverage_chrom.R \\
        --sample      "${meta.id}" \\
        --chrom       "${meta.chrom}" \\
        --tumour-bed  "${tumour_bed}" \\
        --normal-bed  "${normal_bed}" \\
        --segments    "${segment_regions}" \\
        --blacklist   "${blacklist}" \\
        --repeat-bed  "${repeat_bed}" \\
        --outfile     "${prefix}_coverage.rds"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
