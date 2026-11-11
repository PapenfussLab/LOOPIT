process DRAW_CYCLES {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta),
          path(detect_cycles_rds),
          path(annotation),
          path(cnv, stageAs: "cnv/*")

    output:
    tuple val(meta), path("graphs/*"), emit: plots
    path "versions.yml",                   emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def cnv_arg = cnv ? "--cnv \"${cnv}\"" : ""
    """
    mkdir -p graphs

    draw_cycles.R \\
        --sample     "${meta.id}" \\
        --rds        "${detect_cycles_rds}" \\
        --annotation "${annotation}" \\
        ${cnv_arg} \\
        --outdir     "graphs"

    if [ -z "\$(ls -A graphs/*.png 2>/dev/null)" ]; then
        touch graphs/no_cycles_found.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
