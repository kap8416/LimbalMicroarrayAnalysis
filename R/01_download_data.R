## =============================================================================
## 01_download_data.R
## Retrieve GSE38190 (GPL570) and GSE56421 (GPL6244) from NCBI GEO, curate the
## sample set, and write a tidy sample metadata table.
##
## Output:
##   data/raw/<GSE>_series_matrix.rds   normalized expression as provided by GEO
##   data/raw/<GPL>_annotation.rds      probe -> RefSeq / symbol mapping
##   data/raw/sample_metadata.csv       curated 16-sample design
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(GEOquery)
})

## ---- 1. Download series ----------------------------------------------------

fetch_series <- function(gse) {
  rds <- file.path(DIR_RAW, sprintf("%s_series_matrix.rds", gse))
  if (file.exists(rds)) {
    message("cached: ", basename(rds))
    return(readRDS(rds))
  }
  assert_that(!OFFLINE, sprintf("OFFLINE = TRUE but %s is not cached in data/raw/", gse))
  message("downloading ", gse, " from GEO ...")
  es <- GEOquery::getGEO(gse, destdir = DIR_RAW, GSEMatrix = TRUE)[[1]]
  saveRDS(es, rds)
  es
}

es_list <- setNames(lapply(GEO_SERIES$series, fetch_series), GEO_SERIES$series)

## ---- 2. Platform annotation ------------------------------------------------
## Both platforms are annotated to RefSeq accessions, which provides the common
## identifier space for cross-platform integration. Where a probe maps to
## multiple accessions ("A /// B"), the first is retained.

fetch_annotation <- function(gpl) {
  rds <- file.path(DIR_RAW, sprintf("%s_annotation.rds", gpl))
  if (file.exists(rds)) {
    message("cached: ", basename(rds))
    return(readRDS(rds))
  }
  assert_that(!OFFLINE, sprintf("OFFLINE = TRUE but %s annotation is not cached", gpl))
  message("downloading annotation ", gpl, " ...")
  ann <- GEOquery::Table(GEOquery::getGEO(gpl, destdir = DIR_RAW))
  saveRDS(ann, rds)
  ann
}

ann_list <- setNames(lapply(GEO_SERIES$platform, fetch_annotation), GEO_SERIES$platform)

## ---- 3. Curate the sample set ---------------------------------------------
## The full GEO submissions include non-epithelial samples and one conjunctival
## outlier that are removed here. Assignments are explicit rather than parsed
## from GEO titles so that the design is auditable and version-independent.

## Accessions transcribed from the sample map in the analysis scripts that
## produced the published tables (dea_limbal_v3_FINAL_repro_improved.R and
## corto_mra_rewiring_repro_2.R), where they are annotated "verified from GEO".
##
## The full submissions comprise 20 samples; 4 are excluded here, leaving 16.
sample_metadata_full <- tibble::tribble(
  ~gsm,          ~series,     ~platform,  ~group,
  ## ---- GSE56421 / GPL6244 (9 samples) ------------------------------------
  "GSM1361185",  "GSE56421",  "GPL6244",  "Cultured",
  "GSM1361186",  "GSE56421",  "GPL6244",  "Cultured",
  "GSM1361187",  "GSE56421",  "GPL6244",  "Cultured",
  "GSM1361188",  "GSE56421",  "GPL6244",  "Cornea",
  "GSM1361189",  "GSE56421",  "GPL6244",  "Cornea",
  "GSM1361190",  "GSE56421",  "GPL6244",  "Cornea",
  "GSM1361191",  "GSE56421",  "GPL6244",  "Fibroblast",
  "GSM1361192",  "GSE56421",  "GPL6244",  "Fibroblast",
  "GSM1361193",  "GSE56421",  "GPL6244",  "Fibroblast",
  ## ---- GSE38190 / GPL570 (11 samples) ------------------------------------
  "GSM724093",   "GSE38190",  "GPL570",   "Conjunctiva",
  "GSM724094",   "GSE38190",  "GPL570",   "Conjunctiva_outlier",
  "GSM724095",   "GSE38190",  "GPL570",   "Conjunctiva",
  "GSM724096",   "GSE38190",  "GPL570",   "Cornea",
  "GSM724097",   "GSE38190",  "GPL570",   "Cornea",
  "GSM724098",   "GSE38190",  "GPL570",   "Cornea",
  "GSM932369",   "GSE38190",  "GPL570",   "Conjunctiva",
  "GSM932370",   "GSE38190",  "GPL570",   "Limbus",
  "GSM932371",   "GSE38190",  "GPL570",   "Limbus",
  "GSM932372",   "GSE38190",  "GPL570",   "Limbus",
  "GSM932373",   "GSE38190",  "GPL570",   "Limbus"
)

write_result(sample_metadata_full, file.path(DIR_RAW, "sample_metadata_full.csv"))

## Analytical set: fibroblasts and the conjunctival outlier removed.
##
## IMPORTANT — this exclusion must be applied consistently across ALL analytical
## layers. In the original analysis the differential expression step used the
## 16-sample set, but the corto master-regulator and rewiring step was run on a
## 17-sample matrix that retained GSM724094 (see docs/AUDIT.md 3.8). Keeping the
## curated set in one place, here, prevents that divergence from recurring.
sample_metadata <- sample_metadata_full |>
  dplyr::filter(!group %in% c("Fibroblast", "Conjunctiva_outlier")) |>
  dplyr::mutate(
    group     = factor(group, levels = GROUP_LEVELS),
    is_anchor = group == ANCHOR_GROUP
  )

## Verify the accessions against the downloaded GEO records rather than trusting
## the transcription above.
gsm_in_geo <- unlist(lapply(es_list, colnames), use.names = FALSE)
missing_in_geo <- setdiff(sample_metadata_full$gsm, gsm_in_geo)
assert_that(
  length(missing_in_geo) == 0L,
  paste0("These curated GSMs are absent from the downloaded series: ",
         paste(missing_in_geo, collapse = ", "),
         ". Reconcile sample_metadata_full against GEO before proceeding.")
)

if (interactive()) {
  print(lapply(es_list, function(e) data.frame(gsm = colnames(e), title = e$title)))
}

assert_that(
  !any(sample_metadata$gsm %in% EXCLUDED_SAMPLES$gsm),
  "Curated design still contains an excluded sample; check EXCLUDED_SAMPLES."
)
assert_that(
  setequal(setdiff(sample_metadata_full$gsm, sample_metadata$gsm),
           EXCLUDED_SAMPLES$gsm),
  "The samples dropped here do not match EXCLUDED_SAMPLES in R/00_config.R."
)

obs_n <- table(sample_metadata$group)
assert_that(
  all(obs_n[names(EXPECTED_N)] == EXPECTED_N),
  paste0("Group sizes do not match EXPECTED_N. Observed: ",
         paste(names(obs_n), obs_n, sep = "=", collapse = ", "))
)
assert_that(nrow(sample_metadata) == 16L, "Expected 16 curated samples.")

write_result(sample_metadata, file.path(DIR_RAW, "sample_metadata.csv"))
write_result(EXCLUDED_SAMPLES, file.path(DIR_RAW, "excluded_samples.csv"))

message("01_download_data.R complete: ", nrow(sample_metadata), " samples retained.")
