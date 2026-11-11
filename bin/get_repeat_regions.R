#!/usr/bin/env Rscript
# get_repeat_regions.R
# Processes repeat_mask.tsv into a BED file for use in coverage computation
# Called by GET_REPEAT_REGIONS Nextflow process

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option("--reference-dir", type="character", help="Reference directory"),
    make_option("--outfile",       type="character", help="Output repeat mask BED file")
)


opt = parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Processing repeat mask from: ", opt[["reference-dir"]])

config = list(
    output_dir    = dirname(opt$outfile),
    manifest_path = file.path(dirname(opt$outfile), "manifest.csv")
)

result = get_repeat_regions(
    reference_directory  = opt[["reference-dir"]],
    output_directory     = dirname(opt$outfile)
)

if (is.null(result) || !file.exists(result)) {
    stop("get_repeat_regions() did not produce an output file")
}

message("Repeat mask saved to ", result)


