process DETECT_CYCLES {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(coverage), path(breakpoints), path(cnv, stageAs: "cnv/*")

    output:
    tuple val(meta),
          path("${meta.id}_cycles.txt"),
          path("${meta.id}_decomposed_cycles.csv", optional: true),
          path("${meta.id}_filtered_graph.rds"), 
          path("${meta.id}_debug.log"), emit: results
    path "versions.yml",                          emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def purity  = meta.purity ?: 1.0
    def ploidy  = meta.ploidy ?: 2.0
    def mode    = meta.mode   ?: 'SR'
    def cnv_arg = cnv         ? "--cnv \"${cnv}\"" : ""
    def sex_arg = meta.sex    ? "--sex ${meta.sex}" : ""
    """
    detect_cycles.R \\
        --sample            "${meta.id}" \\
        --coverage          "${coverage}" \\
        --breakpoints       "${breakpoints}" \\
        ${cnv_arg} \\
        ${sex_arg} \\
        --mode              "${mode}" \\
        --genome ${params.genome} \\
        --purity            ${purity} \\
        --ploidy            ${ploidy} \\
        --min-cn            "${params.min_cn}" \\
        --min-width         ${params.min_width} \\
        --max-width         ${params.max_width} \\
        --chromosomes       "${params.chromosomes}" \\
        --min-qual-vaf      "${params.min_qual_vaf}" \\
        --min-qual-bp       "${params.min_qual_bp}" \\
        --min-sv-len        "${params.min_sv_len}" \\
        --ecdna-threshold   "${params.ecdna_threshold}" \\
        --outdir            "."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
