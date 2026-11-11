#!/usr/bin/env Rscript
# draw_cycles.R
# Loads CNV (optional), identifies ecDNA cycles, and draws them
# Called by DRAW_CYCLES Nextflow process

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option("--sample",     type="character", help="Sample name"),
    make_option("--rds",        type="character", help="detect_cycles RDS (filtered_edges, merged_nodes, cycles_df, dc, ecdna_threshold, min_cn, cnv_file)"),
    make_option("--annotation", type="character", help="Annotation file"),
    make_option("--outdir",     type="character", default="."),
    make_option("--cnv", type="character", default=NA, help="External CNV file (e.g. FACETS regions)")
)

opt = parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Loading detect_cycles output for ", opt$sample)
filtered_graph = readRDS(opt$rds)

cnv_file = filtered_graph$cnv_file

if (!is.null(filtered_graph$cycles_df) && nrow(filtered_graph$cycles_df) > 0L) {

  #for plotting CNV from an external tool
  if (!is.na(cnv_file) && !is.null(cnv_file) && file.exists(cnv_file)) {
    cnv_df = data.table::fread(cnv_file)
    data.table::setnames(cnv_df, old = "tcn.em", new = "cn", skip_absent = TRUE)
    #convert chrom to character with chr prefix
    cnv_df[, chrom := data.table::fcase(
      as.character(chrom) %in% c("23", "X"), "chrX",
      as.character(chrom) %in% c("24", "Y"), "chrY",
      !startsWith(as.character(chrom), "chr"), paste0("chr", chrom),
      rep(TRUE, .N), as.character(chrom)
    )]
    message("Loaded CNV file: ", cnv_file, " (", nrow(cnv_df), " rows)")
    writeLines(paste("CNV file:", cnv_file, "rows:", nrow(cnv_df)), 
           file.path(opt$outdir, "cnv_debug.txt"))
  } else {
    cnv_df = data.frame()
  }

  ecdna_cycle_ids = as.character(
    filtered_graph$dc[filtered_graph$dc$ecdna == TRUE, "cycle_id"][[1]]
  )

  ecdna_cycles_df = filtered_graph$cycles_df[
    filtered_graph$cycles_df$cycle_id %in% ecdna_cycle_ids,
  ]

  message("ecDNA cycle IDs: ", paste(ecdna_cycle_ids, collapse=", "))
  message("ecdna_cycles_df rows: ", nrow(ecdna_cycles_df))

  draw_cycles(filtered_graph$filtered_edges, filtered_graph$merged_nodes,
              ecdna_cycles_df, opt$annotation, cnv_df,
              opt$outdir, filtered_graph$min_cn)

} else {
  message("No cycles found — skipping drawing")
}