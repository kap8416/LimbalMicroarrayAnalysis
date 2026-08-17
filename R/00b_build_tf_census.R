## =============================================================================
## 00b_build_tf_census.R
## Build the transcription-factor census used as corto centroids, and write it to
## data/raw/human_tfs.txt.
##
## WHY THIS EXISTS
## Master regulator analysis works by inferring a regulon for every candidate
## regulator and then testing which regulons are coordinately shifted between
## conditions. The candidate set therefore determines what the analysis can find,
## and it must be defined independently of the biology under study.
##
## The original analysis used a hand-written list of 112 "TF candidates", assembled
## as "known limbal-relevant TFs" plus a partial catalogue (see
## provenance/corto_mra_rewiring_repro_2.R, line 241). Two problems follow
## (docs/AUDIT.md 3.9):
##
##   1. CIRCULARITY. Genes were selected for their known relevance to the limbus,
##      then reported as the regulators the analysis discovered.
##
##   2. NON-TRANSCRIPTION-FACTORS. The list includes secreted ligands (BMP4, SHH,
##      WNT5A, TGFB1, DLL1, JAG1), receptors (FZD7, PTCH1), kinases (CDK4, IKBKB),
##      an intermediate filament (VIM), a Ras GTPase-activating protein (NF1), a
##      scaffold protein (NF2), RNA-binding proteins (YBX1, HNRNPK, SRSF1), an E3
##      ligase (MDM2), and "AP1", which is not an HGNC symbol. Several of these are
##      reported in the manuscript as master transcriptional regulators.
##
## This script builds a census from an established resource instead. Preference
## order:
##   1. data/raw/TF_names_v_1.01.txt   — Lambert et al. (2018) human TF census,
##      the standard reference (~1,639 factors). Download once:
##        curl -L -o data/raw/TF_names_v_1.01.txt \
##          http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
##   2. Gene Ontology molecular-function terms via org.Hs.eg.db, which needs no
##      download but is broader and noisier.
##
## Output:
##   data/raw/human_tfs.txt          one symbol per line
##   data/raw/human_tfs_source.txt   which resource was used, and the count
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

LAMBERT_FILE <- file.path(DIR_RAW, "TF_names_v_1.01.txt")
OUT_FILE     <- file.path(DIR_RAW, "human_tfs.txt")

## GO molecular-function terms describing transcriptional regulation.
## GO:0003700  DNA-binding transcription factor activity
## GO:0000981  DNA-binding transcription factor activity, RNA polymerase II-specific
## GO:0140110  transcription regulator activity
GO_TF_TERMS <- c("GO:0003700", "GO:0000981", "GO:0140110")

if (file.exists(LAMBERT_FILE)) {
  tfs <- unique(trimws(readLines(LAMBERT_FILE, warn = FALSE)))
  tfs <- tfs[nzchar(tfs) & !grepl("^#", tfs)]
  src <- sprintf("Lambert et al. 2018 human TF census (%s)", basename(LAMBERT_FILE))
} else {
  message("Lambert census not found at ", LAMBERT_FILE)
  message("Falling back to Gene Ontology molecular-function terms.")
  message("For a publication-grade run, download the census:")
  message("  curl -L -o ", LAMBERT_FILE,
          " http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt")

  eg <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = GO_TF_TERMS, keytype = "GOALL", columns = "SYMBOL"
  ))
  tfs <- sort(unique(eg$SYMBOL[!is.na(eg$SYMBOL)]))
  src <- sprintf("Gene Ontology GOALL {%s} via org.Hs.eg.db",
                 paste(GO_TF_TERMS, collapse = ", "))
}

message(sprintf("\ncandidate regulators: %d", length(tfs)))
message("source: ", src)

## A census far below ~1,000 factors indicates a truncated or hand-picked list,
## which is what produced the circularity described above. Fail rather than let a
## small set silently shape the master regulator results.
assert_that(
  length(tfs) >= 800L,
  sprintf(paste0(
    "Only %d candidate regulators were assembled; a human TF census holds roughly ",
    "1,600. A set this small biases master regulator analysis toward whatever it ",
    "happens to contain. Download the Lambert census into data/raw/ and re-run:\n",
    "  curl -L -o %s http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"),
    length(tfs), LAMBERT_FILE)
)

## Report which regulators named in the manuscript are, and are not, in the census.
MANUSCRIPT_REGULATORS <- c(
  "ETS1", "THRB", "ZEB1", "ZEB2", "PPARG", "MEF2C", "MITF", "IKZF1", "FLI1",
  "NFIX", "SMAD5", "KLF4", "MAFB", "EGR1", "PITX2", "RB1", "PAX6", "OVOL1",
  "PTCH1", "HMGA2", "ELF3",
  ## Named in the manuscript as regulators but not transcription factors:
  "BMP4", "FZD7", "VIM", "NF1"
)
present <- MANUSCRIPT_REGULATORS %in% tfs
if (any(!present)) {
  message("\nRegulators named in the manuscript that are ABSENT from the TF census:")
  message("  ", paste(MANUSCRIPT_REGULATORS[!present], collapse = ", "))
  message("  These are not transcription factors and cannot be master ",
          "transcriptional regulators. Any NES reported for them came from ",
          "including them in a hand-written centroid list.")
}

writeLines(tfs, OUT_FILE)
writeLines(
  c(sprintf("source: %s", src),
    sprintf("n_factors: %d", length(tfs)),
    sprintf("built: %s", as.character(Sys.time())),
    sprintf("absent_manuscript_regulators: %s",
            paste(MANUSCRIPT_REGULATORS[!present], collapse = ", "))),
  file.path(DIR_RAW, "human_tfs_source.txt")
)

message("\nwrote ", OUT_FILE, " (", length(tfs), " factors)")
message("00b_build_tf_census.R complete.")
