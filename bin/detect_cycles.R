#!/usr/bin/env Rscript
# detect_cycles.R
# Runs SCC detection, cycle decomposition and scoring
# Called by DETECT_CYCLES Nextflow process
# Accepts pre-computed coverage, breakpoints and segments

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option("--sample",          type="character", help="Sample name"),
    make_option("--coverage",        type="character", help="Merged coverage RDS (list tumour+normal)"),
    make_option("--breakpoints",     type="character", help="Breakpoints RDS"),
    make_option("--cnv",             type="character", default=NA,    help="CNV file (optional)"),
    make_option("--mode",            type="character", default="SR",  help="SR or LR"),
    make_option("--genome",          type="character", default="GRCh38"),
    make_option("--purity",          type="double",    default=1.0,   help="Tumour purity"),
    make_option("--ploidy",          type="double",    default=2.0,   help="Tumour ploidy"),
    make_option("--min-cn",          type="character", default="2,4"),
    make_option("--min-width",       type="integer",   default=10000),
    make_option("--max-width",       type="integer",   default=1e7),
    make_option("--chromosomes",     type="character", default="NA"),
    make_option("--min-qual-vaf",    type="character", default="AUTO"),
    make_option("--min-qual-bp",     type="character", default="50"),
    make_option("--min-sv-len",      type="character", default="300"),
    make_option("--ecdna-threshold", type="character", default="2.5"),
    make_option("--outdir",          type="character", default="."),
    make_option("--sex", type="character", default=NA, help="Sample sex (M/F). If not provided, inferred from chrY/autosome coverage ratio.")
)

opt = parse_args(OptionParser(option_list = option_list))


mode_defaults = list(
  SR = list(
    min_qual_vaf = "0.2,0.4,0.6"
  ),
  LR = list(
    min_qual_vaf = "0.3,0.5,0.6,0.8"
  )
)


defaults_for_mode = mode_defaults[[opt$mode]]
if (is.null(defaults_for_mode)) {
  message("Unrecognised mode '", opt$mode, "' — falling back to SR filter defaults")
  defaults_for_mode = mode_defaults[["SR"]]
}


min_qual_vaf_val = if (!is.null(opt[["min-qual-vaf"]]) && opt[["min-qual-vaf"]] != "AUTO") {
  opt[["min-qual-vaf"]]
} else {
  defaults_for_mode$min_qual_vaf
}

debug_log = file.path(opt$outdir, paste0(opt$sample, "_debug.log"))
log_con   = file(debug_log, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Running cycle detection for ", opt$sample)

#load pre-computed inputs
message("Loading coverage...")
coverage = readRDS(opt$coverage)

message("Loading breakpoints...")
breakpoints = readRDS(opt$breakpoints)

cnv_file    = if (is.na(opt$cnv) || opt$cnv == "NA" || opt$cnv == "") NA else opt$cnv
chromosomes = if (opt$chromosomes == "NA") NA else strsplit(opt$chromosomes, ",")[[1]]

parse_values = function(x) as.numeric(strsplit(x, ",")[[1]])

filter_grid = expand.grid(
    MIN_QUAL_BREAKPOINT = parse_values(opt[["min-qual-bp"]]),
    MIN_QUAL_VAF        = parse_values(min_qual_vaf_val),
    MIN_SV_LEN          = parse_values(opt[["min-sv-len"]]),
    MIN_CN              = parse_values(opt[["min-cn"]]),
    ECDNA_THRESHOLD     = parse_values(opt[["ecdna-threshold"]])
)

#calculate cn
all_cn = calculate_cn(coverage$normal, coverage$tumour)

sex = infer_sample_sex(coverage$tumour)

if (!is.na(opt$sex) && opt$sex %in% c("M", "F") && opt$sex != sex) {
  message("NOTE: clinical sex annotation (", opt$sex,
          ") differs from coverage-inferred sex (", sex,
          ") using coverage-inferred sex for CN baseline; ",
          "this discordance may indicate LOY or a sample-labelling issue worth reviewing")
} else if (!is.na(opt$sex) && opt$sex %in% c("M", "F")) {
  message("Coverage-inferred sex (", sex, ") agrees with clinical annotation")
}


if (!is.na(opt$purity) && !is.na(opt$ploidy)) {
  message("Rescaling copy number")
all_cn = rescale_cn(all_cn, opt$purity, opt$ploidy, opt$genome, sex)
}

#run cycle detection

result = detect_cycles(
    outdir          = opt$outdir,
    all_cn          = all_cn,
    breakpoints_df  = breakpoints$breakpoint_df,
    min_width       = opt[["min-width"]],
    chromosomes     = chromosomes,
    filter_grid     = filter_grid,
    mode            = opt$mode,
    ploidy          = opt$ploidy
)




cycles_file = file.path(opt$outdir, paste0(opt$sample, "_cycles.txt"))
writeLines(result$cycles_out, cycles_file)

message("Cycles written to: ", cycles_file)

dc_file = file.path(opt$outdir, paste0(opt$sample, "_decomposed_cycles.csv"))

if (!is.null(result$cycles_df) && nrow(result$cycles_df) > 0) {

  if (!is.null(result$dc) && nrow(result$dc) > 0) {
    dc = result$dc
  } else {
    dc = decompose_cycles(result$cycles_df,
                          merged_nodes = result$merged_nodes,
                          threshold    = result$ecdna_threshold)
  }

  if (!is.null(dc) && nrow(dc) > 0) {
    data.table::fwrite(dc, dc_file, sep = ",", col.names = TRUE)
    message("Decomposed cycles written to: ", dc_file)
  } else {
    message("No decomposed cycles to write")
    file.create(dc_file)
  }

} else {
  message("No cycles found — skipping decomposition")
  dc = data.table::data.table()
  file.create(dc_file)
}


saveRDS(
  list(
    filtered_edges  = result$filtered_edges,
    merged_nodes    = result$merged_nodes,
    cycles_df       = result$cycles_df,
    dc              = dc,
    ecdna_threshold = result$ecdna_threshold,
    min_cn          = result$min_cn,
    cnv_file        = cnv_file    
  ),
  file.path(opt$outdir, paste0(opt$sample, "_filtered_graph.rds"))
)

sink(type = "message")
sink(type = "output")
close(log_con)


