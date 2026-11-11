#!/usr/bin/env Rscript
# get_breakpoints.R
# Parses SV VCF and filters breakpoints using repeat mask
# Called by GET_BREAKPOINTS Nextflow process

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option("--sample",     type="character", help="Sample name"),
    make_option("--vcf",        type="character", help="SV VCF file (GRIDSS or Sniffles)"),
    make_option("--mode",       type="character", default="SR", help="SR or LR"),
    make_option("--repeat-bed", type="character", help="Repeat mask BED file"),
    make_option("--outdir",     type="character", default="."),
    make_option("--outfile",    type="character", help="Output breakpoints RDS file")
)

opt = parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Getting breakpoints for ", opt$sample)
message("Mode: ", opt$mode)
message("VCF: ", opt$vcf)

config = list(
    output_dir    = opt$outdir,
    manifest_path = file.path(opt$outdir, "manifest.csv")
)

breakpoints = get_breakpoints(
    sample_name      = opt$sample,
    vcf_file         = opt$vcf,
    mode             = opt$mode,
    repeat_file      = opt[["repeat-bed"]],
    output_directory = opt$outdir,
    chromosomes      = NA
)

saveRDS(breakpoints, opt$outfile)
message("Breakpoints saved to ", opt$outfile)
message("  breakpoint_df rows: ", nrow(breakpoints$breakpoint_df))
message("  breakend_df rows:   ", nrow(breakpoints$breakend_df))
