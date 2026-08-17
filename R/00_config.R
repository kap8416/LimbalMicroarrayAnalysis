## =============================================================================
## 00_config.R
## Global configuration: paths, sample design, thresholds, palettes.
## Sourced by every downstream script. Modify here, not in the analysis scripts.
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---- Project paths ---------------------------------------------------------
## Assumes the working directory is the repository root.
PROJ        <- normalizePath(".", mustWork = TRUE)
DIR_RAW     <- file.path(PROJ, "data", "raw")
DIR_PROC    <- file.path(PROJ, "data", "processed")
DIR_RES     <- file.path(PROJ, "results")
DIR_DE      <- file.path(DIR_RES, "de")
DIR_ENR     <- file.path(DIR_RES, "enrichment")
DIR_PPI     <- file.path(DIR_RES, "ppi")
DIR_MRA     <- file.path(DIR_RES, "mra")
DIR_FIG     <- file.path(PROJ, "figures", "output")

for (d in c(DIR_RAW, DIR_PROC, DIR_DE, DIR_ENR, DIR_PPI, DIR_MRA, DIR_FIG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---- Reproducibility -------------------------------------------------------
SEED <- 42
set.seed(SEED)

## Set to TRUE to reuse cached GEO downloads in data/raw/ instead of
## re-querying GEO (useful on machines without network access).
OFFLINE <- FALSE

## ---- Datasets --------------------------------------------------------------
GEO_SERIES <- tibble::tribble(
  ~series,     ~platform,  ~platform_chip,
  "GSE38190",  "GPL570",   "Affymetrix HG-U133 Plus 2.0",
  "GSE56421",  "GPL6244",  "Affymetrix HuGene-1_0-st"
)

## Samples excluded during curation, with justification.
##  - limbal fibroblasts: non-epithelial, outside the scope of the contrasts.
##  - GSM724094: conjunctival outlier; within-group Pearson r = 0.92 versus
##    r = 0.98-0.99 for the remaining conjunctival replicates.
EXCLUDED_SAMPLES <- tibble::tribble(
  ~gsm,          ~reason,
  "GSM1361191",  "limbal fibroblast (non-epithelial)",
  "GSM1361192",  "limbal fibroblast (non-epithelial)",
  "GSM1361193",  "limbal fibroblast (non-epithelial)",
  "GSM724094",   "conjunctival outlier (within-group r = 0.92)"
)

## Analytical groups. Cornea is the only group present on both platforms and
## therefore serves as the anchor for cross-platform harmonization.
GROUP_LEVELS <- c("Limbus", "Cornea", "Conjunctiva", "Cultured")
ANCHOR_GROUP <- "Cornea"

## Expected final design (asserted in 02_integration.R).
EXPECTED_N <- c(Limbus = 4L, Cornea = 6L, Conjunctiva = 3L, Cultured = 3L)

## ---- Integration QC --------------------------------------------------------
N_TRANSCRIPTS_EXPECTED <- 45378L   # RefSeq-level estimates before collapse
N_GENES_EXPECTED       <- 24002L   # gene-level estimates after MAD collapse
CONCORDANCE_MIN        <- 0.95     # pre-specified cornea cross-platform r

## ---- Contrast definitions --------------------------------------------------
## Limbus is the reference condition in every contrast, so a positive log2FC
## means higher expression in limbus.
##
## Fold-change thresholds are contrast-specific because platform confounding
## differs: C2 compares samples from a single platform and tolerates a lower
## threshold, whereas C1 and C3 are cross-platform.
##
## C3 is FULLY CONFOUNDED with platform (limbus = GPL570, cultured = GPL6244)
## and is reported as exploratory throughout.
CONTRASTS <- tibble::tribble(
  ~id,   ~reference, ~comparator,    ~lfc_cut, ~fdr_cut, ~cross_platform, ~exploratory, ~label,
  "C1",  "Limbus",   "Cornea",       1.0,      0.05,     TRUE,            FALSE,        "Limbus vs cornea",
  "C2",  "Limbus",   "Conjunctiva",  0.5,      0.05,     FALSE,           FALSE,        "Limbus vs conjunctiva",
  "C3",  "Limbus",   "Cultured",     1.0,      0.05,     TRUE,            TRUE,         "Limbus vs cultured limbal epithelial cells"
)

## ---- PPI parameters --------------------------------------------------------
STRING_VERSION   <- "12.0"
STRING_API       <- "https://string-db.org/api"
STRING_SPECIES   <- 9606L
STRING_SCORE_MIN <- 700L          # combined confidence >= 0.70, STRING scale 0-1000
STRING_NETWORK   <- "functional"
N_HUBS_REPORT    <- 15L

## Interaction source: "local" or "api".
##
## "local" is the default and the recommended setting. It reads the STRING flat
## files, which contain the complete interactome with no query-size limit, and is
## exactly reproducible because the files are versioned.
##
## "api" queries the REST endpoint in a SINGLE request. It is only valid when every
## contrast's gene set fits within STRING_API_MAX identifiers; 05_ppi_networks.R
## stops rather than chunking. Chunking is prohibited because the /network endpoint
## returns only interactions among the submitted identifiers, so splitting the list
## silently discards every interaction spanning two chunks — the defect that left
## the original C1 and C3 networks missing 12% and 5% of their edges
## (docs/AUDIT.md 3.3).
STRING_SOURCE        <- "local"
STRING_API_MAX       <- 1900L
STRING_LOCAL_LINKS   <- "9606.protein.links.v12.0.txt.gz"
STRING_LOCAL_ALIASES <- "9606.protein.aliases.v12.0.txt.gz"

## Maximum number of maximal cliques to enumerate for exact MCC. Above this the
## metric is dropped from the composite hub score and the omission is recorded.
## A degree-based surrogate is never substituted: degree^2 is a monotone function
## of degree and would double-weight a term already in the composite
## (docs/AUDIT.md 3.4).
MCC_CLIQUE_BUDGET <- 5e6

## ---- MRA / co-expression parameters ---------------------------------------
CORTO_NTOP       <- 12000L         # top most-variable genes used as input
CORTO_NBOOTSTRAP <- 1000L
CORTO_PVALUE     <- 1e-6
MRA_MINSIZE      <- 10L           # minimum regulon size
REWIRE_DELTA_RHO <- 0.6           # |delta Spearman rho| threshold for a rewired edge

## ---- Enrichment parameters -------------------------------------------------
GO_ONTOLOGY  <- "BP"
KEGG_RELEASE <- "117.0"
ENR_FDR_CUT  <- 0.05
ENR_ORGDB    <- "org.Hs.eg.db"

## ---- Figure palette (shared with figures/scripts/fig_style.py) ------------
PAL_GROUP <- c(
  Limbus      = "#1B7A6E",   # teal
  Cornea      = "#C4622D",   # burnt orange
  Conjunctiva = "#3D5A8A",   # slate blue
  Cultured    = "#8A6BA8"    # muted violet
)
PAL_DIRECTION <- c(Up = "#B03A2E", Down = "#2471A3", NS = "#BFBFBF")
PAL_CONTRAST  <- c(C1 = "#C4622D", C2 = "#3D5A8A", C3 = "#8A6BA8")

## ---- Helpers ---------------------------------------------------------------

#' Write a results table with consistent formatting.
#'
#' Any existing file is removed first. Outputs left over from an earlier run can
#' carry restrictive permissions, macOS immutable/locked flags, or an open handle
#' from a spreadsheet application, and `readr::write_csv` then fails with only
#' "Cannot open file for writing", which does not say why. Unlinking first handles
#' the common cases, and the fallback below reports the actual cause.
write_result <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(path)) {
    suppressWarnings(try(Sys.chmod(path, "0644"), silent = TRUE))
    if (!file.remove(path) && file.exists(path)) {
      stop(sprintf(paste0(
        "Cannot replace an existing output file:\n  %s\n",
        "It is present but could not be removed. Usual causes:\n",
        "  - the file is open in Excel, Numbers or another application: close it;\n",
        "  - it is locked or flagged immutable on macOS:\n",
        "        chmod u+w '%s' ; chflags nouchg '%s'\n",
        "  - the directory is not writable by this user.\n",
        "Stale outputs from a previous run can also simply be deleted:\n",
        "        rm -f results/de/*.csv results/ppi/*.csv results/mra/*.csv"),
        path, path, path), call. = FALSE)
    }
  }

  ok <- tryCatch({ readr::write_csv(x, path); TRUE },
                 error = function(e) { message("  write failed: ", conditionMessage(e)); FALSE })
  assert_that(ok, sprintf("Could not write %s", path))

  message(sprintf("  wrote %s (%d rows)", basename(path), nrow(x)))
  invisible(path)
}

#' Hard assertion with an informative message; used for QC gates.
assert_that <- function(condition, msg) {
  if (!isTRUE(condition)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

#' Median absolute deviation across a numeric vector, NA-safe.
row_mad <- function(m) apply(m, 1, stats::mad, na.rm = TRUE)

message("00_config.R loaded | seed = ", SEED, " | offline = ", OFFLINE)
