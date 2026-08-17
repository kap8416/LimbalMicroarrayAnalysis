## =============================================================================
## 99_ingest_legacy_results.R
## Utility: map result files produced by the original exploratory analysis
## (TranscriptomadeOjo.R) onto the canonical file names expected by the figure
## scripts, so figures can be regenerated without re-running the full pipeline.
##
## Renaming performed:
##   centrality_Limbus_vs_Cornea.csv       -> results/ppi/centrality_C1.csv
##   centrality_Limbus_vs_Conjunctiva.csv  -> results/ppi/centrality_C2.csv
##   centrality_Limbus_vs_AMLE.csv         -> results/ppi/centrality_C3.csv
##   DEA_Limbus_vs_AMLE.csv                -> results/de/DEA_C3_Limbus_vs_Cultured.csv
##   bipartite_*_v3.module_membership.csv  -> results/de/deg_membership.csv
##
## "AMLE" in the legacy filenames denotes the cultured limbal epithelial cell
## group, referred to as "Cultured" throughout this repository and as
## "cultured limbal epithelial cells" in the manuscript.
##
## Usage:  Rscript R/99_ingest_legacy_results.R <legacy_dir>
## =============================================================================

source("R/00_config.R")

args   <- commandArgs(trailingOnly = TRUE)
legacy <- if (length(args) >= 1L) args[[1]] else file.path(DIR_RAW, "legacy")
assert_that(dir.exists(legacy), sprintf("Legacy directory not found: %s", legacy))

copy_if <- function(src, dest) {
  src <- file.path(legacy, src)
  if (!file.exists(src)) {
    message("  missing (skipped): ", basename(src))
    return(invisible(FALSE))
  }
  file.copy(src, dest, overwrite = TRUE)
  message("  ", basename(src), " -> ", basename(dest))
  invisible(TRUE)
}

message("ingesting PPI centrality tables ...")
copy_if("centrality_Limbus_vs_Cornea.csv",      file.path(DIR_PPI, "centrality_C1.csv"))
copy_if("centrality_Limbus_vs_Conjunctiva.csv", file.path(DIR_PPI, "centrality_C2.csv"))
copy_if("centrality_Limbus_vs_AMLE.csv",        file.path(DIR_PPI, "centrality_C3.csv"))

message("ingesting differential expression tables ...")
copy_if("DEA_Limbus_vs_AMLE.csv",     file.path(DIR_DE, "DEA_C3_Limbus_vs_Cultured.csv"))
copy_if("DEA_Limbus_vs_Cornea.csv",   file.path(DIR_DE, "DEA_C1_Limbus_vs_Cornea.csv"))
copy_if("DEA_Limbus_vs_Conjunctiva.csv", file.path(DIR_DE, "DEA_C2_Limbus_vs_Conjunctiva.csv"))

## ---- Cross-contrast membership --------------------------------------------
## The legacy bipartite table already carries per-contrast status, logFC and FDR;
## only the column names differ from the canonical schema.

mm_src <- file.path(legacy, "bipartite_network_clean_v3.module_membership.csv")
if (file.exists(mm_src)) {
  mm <- readr::read_csv(mm_src, show_col_types = FALSE) |>
    dplyr::rename(
      Symbol      = gene,
      C1_status   = C1, C1_logFC = C1_logFC, C1_adj.P.Val = C1_FDR,
      C2_status   = C2, C2_logFC = C2_logFC, C2_adj.P.Val = C2_FDR,
      C3_status   = C3, C3_logFC = C3_logFC, C3_adj.P.Val = C3_FDR
    ) |>
    ## Legacy status labels are upper case.
    dplyr::mutate(dplyr::across(dplyr::ends_with("_status"),
                                ~ dplyr::recode(.x, UP = "Up", DOWN = "Down", NS = "NS"))) |>
    dplyr::mutate(
      n_up = (C1_status == "Up") + (C2_status == "Up") + (C3_status == "Up")
    )
  write_result(mm, file.path(DIR_DE, "deg_membership.csv"))
} else {
  message("  missing (skipped): bipartite_network_clean_v3.module_membership.csv")
}

## ---- Counts summary --------------------------------------------------------
## Two legacy counts files exist. The transcript-level one (45,378 features,
## 351/29 - 76/39 - 51/6) PRECEDES the isoform collapse and does NOT correspond
## to the published results; only the gene-level file is ingested.

for (f in c("DEA_counts_summary-effc57a2.csv", "DEA_counts_summary.csv")) {
  p <- file.path(legacy, f)
  if (!file.exists(p)) next
  cs <- readr::read_csv(p, show_col_types = FALSE)
  if ("total_raw" %in% names(cs) && max(cs$total_raw, na.rm = TRUE) == 45378) {
    message("  ", f, " is transcript-level (pre-collapse); NOT ingested")
    next
  }
  write_result(cs, file.path(DIR_DE, "DEA_counts_summary.csv"))
  break
}

message("99_ingest_legacy_results.R complete.")
