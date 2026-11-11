#!/usr/bin/env Rscript
# get_segments.R
# Defines segment boundaries from breakpoints and reference regions
# Called by GET_SEGMENTS Nextflow process

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option("--sample",            type="character", help="Sample name"),
    make_option("--breakpoints",       type="character", help="Breakpoints RDS file"),
    make_option("--reference-regions", type="character", help="Reference regions RDS file"),
    make_option("--chromosomes",       type="character", default="NA"),
    make_option("--outfile",           type="character", help="Output segments BED file")
)

opt = parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Getting segments for ", opt$sample)

outdir = dirname(opt$outfile)
config = list(
    output_dir    = outdir,
    manifest_path = file.path(outdir, "manifest.csv")
)

#load pre-computed inputs
breakpoints       = readRDS(opt$breakpoints)
reference_regions = readRDS(opt[["reference-regions"]])

chromosomes = if (opt$chromosomes == "NA") NA else strsplit(opt$chromosomes, ",")[[1]]

result = get_segments(
    region_ct        = reference_regions$region_ct,
    chr_regions      = reference_regions$chr_regions,
    breakpoints      = breakpoints$breakpoint_bed,
    sample_name      = opt$sample,
    output_directory = outdir,
    chromosomes      = chromosomes
)

if (!is.null(result) && file.exists(result) && result != opt$outfile) {
    file.copy(result, opt$outfile, overwrite = TRUE)
}

message("Segments saved to ", opt$outfile)
