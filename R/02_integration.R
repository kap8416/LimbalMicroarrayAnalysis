## =============================================================================
## 02_integration.R
## Cross-platform harmonization of GSE38190 (GPL570) and GSE56421 (GPL6244).
##
## Strategy (three steps, because corneal epithelium is the only group present in
## both datasets and therefore the only available cross-platform anchor):
##   (1) within-platform quantile normalization  [limma::normalizeBetweenArrays]
##   (2) cornea-anchored cross-platform shift correction
##   (3) residual batch correction               [sva::ComBat, mod = ~is_Cornea]
##
## Probes are first summarized at RefSeq-transcript level, then collapsed to
## gene level by retaining, per gene symbol, the transcript with the highest
## median absolute deviation (MAD) across samples.
##
## Output:
##   data/processed/expr_transcript.rds   RefSeq-level matrix  (45,378 x 16)
##   data/processed/expr_gene.rds         gene-level matrix    (24,002 x 16)
##   data/processed/expr_gene.csv         same, for the Python figure scripts
##   results/de/integration_qc.csv        concordance and PCA variance
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(limma)
  library(sva)
})

meta <- readr::read_csv(file.path(DIR_RAW, "sample_metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(group = factor(group, levels = GROUP_LEVELS))

## ---- 1. Probe -> RefSeq summarization --------------------------------------

#' Map probes to RefSeq accessions and average replicate probes per accession.
#'
#' @param expr numeric matrix, probes x samples (log2 scale as provided by GEO)
#' @param ann  platform annotation table with probe ID and RefSeq columns
summarize_to_refseq <- function(expr, ann, probe_col, refseq_col) {
  ann <- ann[!is.na(ann[[refseq_col]]) & ann[[refseq_col]] != "", , drop = FALSE]
  keep <- intersect(rownames(expr), as.character(ann[[probe_col]]))
  expr <- expr[keep, , drop = FALSE]
  ann  <- ann[match(keep, as.character(ann[[probe_col]])), , drop = FALSE]

  ## Multi-mapping probes ("NM_1 /// NM_2") are assigned to their first accession.
  refseq <- vapply(strsplit(as.character(ann[[refseq_col]]), " /// ", fixed = TRUE),
                   `[`, character(1), 1)

  ## Average over probes sharing an accession.
  out <- rowsum(expr, group = refseq, reorder = TRUE, na.rm = TRUE)
  n   <- as.vector(table(refseq)[rownames(out)])
  sweep(out, 1, n, "/")
}

#' Locate a column by any of several candidate names.
#'
#' GEO platform annotation tables are not consistently named across releases
#' (GPL570 uses "RefSeq Transcript ID", GPL6244 has used "GB_ACC", "gene_assignment"
#' and "RefSeq"). Failing loudly here beats silently annotating against the wrong
#' column, which would corrupt every downstream identifier.
pick_col <- function(ann, candidates, what, platform) {
  hit <- intersect(candidates, names(ann))
  assert_that(
    length(hit) > 0L,
    sprintf("%s: no %s column found. Candidates tried: %s. Available: %s",
            platform, what, paste(candidates, collapse = ", "),
            paste(names(ann), collapse = ", "))
  )
  hit[1]
}

load_platform <- function(series, platform, probe_col = NULL, refseq_col = NULL) {
  es  <- readRDS(file.path(DIR_RAW, sprintf("%s_series_matrix.rds", series)))
  ann <- readRDS(file.path(DIR_RAW, sprintf("%s_annotation.rds", platform)))

  if (is.null(probe_col)) {
    probe_col <- pick_col(ann, c("ID", "probeset_id", "PROBEID"), "probe ID", platform)
  }
  if (is.null(refseq_col)) {
    refseq_col <- pick_col(
      ann,
      c("RefSeq Transcript ID", "RefSeq", "GB_ACC", "GB_LIST", "RefSeq_ID"),
      "RefSeq accession", platform
    )
  }
  message(sprintf("  %s: probe column '%s', RefSeq column '%s'",
                  platform, probe_col, refseq_col))

  expr <- Biobase::exprs(es)
  ## Some GEO series matrices are supplied on the linear scale.
  if (max(expr, na.rm = TRUE) > 100) {
    message("  ", series, ": values > 100 detected, applying log2(x + 1)")
    expr <- log2(expr + 1)
  }

  gsm  <- meta$gsm[meta$series == series]
  assert_that(all(gsm %in% colnames(expr)),
              sprintf("%s: curated GSMs missing from the series matrix.", series))
  expr <- expr[, gsm, drop = FALSE]

  ## (1) within-platform quantile normalization
  expr <- limma::normalizeBetweenArrays(expr, method = "quantile")

  summarize_to_refseq(expr, ann, probe_col, refseq_col)
}

message("summarizing GPL570 (GSE38190) ...")
m570 <- load_platform("GSE38190", "GPL570")

message("summarizing GPL6244 (GSE56421) ...")
m6244 <- load_platform("GSE56421", "GPL6244")

## ---- 2. Cornea-anchored cross-platform shift correction --------------------
## Restricted to transcripts measured on both platforms. The corneal group is
## profiled in both datasets, so the median cornea-to-cornea difference per
## transcript estimates the platform offset and is subtracted from GPL6244.

common <- intersect(rownames(m570), rownames(m6244))
message("transcripts shared across platforms: ", length(common))
m570  <- m570[common, , drop = FALSE]
m6244 <- m6244[common, , drop = FALSE]

cornea_570  <- meta$gsm[meta$series == "GSE38190" & meta$group == "Cornea"]
cornea_6244 <- meta$gsm[meta$series == "GSE56421" & meta$group == "Cornea"]

shift <- rowMeans(m570[, cornea_570, drop = FALSE], na.rm = TRUE) -
         rowMeans(m6244[, cornea_6244, drop = FALSE], na.rm = TRUE)
m6244_shifted <- sweep(m6244, 1, shift, "+")

expr_tx <- cbind(m570, m6244_shifted)[, meta$gsm, drop = FALSE]

## ---- 3. Residual ComBat correction ----------------------------------------
## The biological covariate protected during batch correction is corneal
## identity, the only group observed in both batches. Protecting the full group
## factor is not possible because group and batch are partially confounded.

batch     <- factor(meta$series)
is_Cornea <- as.integer(meta$group == "Cornea")
mod       <- model.matrix(~ is_Cornea)

expr_tx <- sva::ComBat(dat = as.matrix(expr_tx), batch = batch,
                       mod = mod, par.prior = TRUE, prior.plots = FALSE)

if (nrow(expr_tx) != N_TRANSCRIPTS_EXPECTED) {
  warning(sprintf(
    "Transcript count is %d; the manuscript reports %d. Check platform annotation versions.",
    nrow(expr_tx), N_TRANSCRIPTS_EXPECTED))
}

saveRDS(expr_tx, file.path(DIR_PROC, "expr_transcript.rds"))

## ---- 4. Transcript -> gene collapse by maximum MAD -------------------------
## Retaining the most variable isoform per gene preserves the transcript that
## carries the informative signal, rather than diluting it by averaging isoforms
## with divergent expression behaviour.

ann570  <- readRDS(file.path(DIR_RAW, "GPL570_annotation.rds"))
rs_col  <- pick_col(ann570, c("RefSeq Transcript ID", "RefSeq", "GB_ACC"),
                    "RefSeq accession", "GPL570")
sym_col <- pick_col(ann570, c("Gene Symbol", "GENE_SYMBOL", "Symbol"),
                    "gene symbol", "GPL570")

tx2sym <- ann570 |>
  dplyr::transmute(
    refseq = vapply(strsplit(as.character(.data[[rs_col]]),  " /// ", fixed = TRUE),
                    `[`, character(1), 1),
    symbol = vapply(strsplit(as.character(.data[[sym_col]]), " /// ", fixed = TRUE),
                    `[`, character(1), 1)
  ) |>
  dplyr::filter(!is.na(refseq), refseq != "", !is.na(symbol), symbol != "") |>
  dplyr::distinct(refseq, .keep_all = TRUE)

## Transcripts without a symbol keep their accession as identifier, which is why
## predicted XM_/XR_ accessions persist in the gene-level matrix.
sym <- tx2sym$symbol[match(rownames(expr_tx), tx2sym$refseq)]
sym[is.na(sym)] <- rownames(expr_tx)[is.na(sym)]

mads <- row_mad(expr_tx)
keep <- tapply(seq_len(nrow(expr_tx)), sym, function(i) i[which.max(mads[i])])
expr_gene <- expr_tx[unlist(keep), , drop = FALSE]
rownames(expr_gene) <- names(keep)

if (nrow(expr_gene) != N_GENES_EXPECTED) {
  warning(sprintf(
    "Gene count is %d; the manuscript reports %d. Check annotation versions.",
    nrow(expr_gene), N_GENES_EXPECTED))
}

saveRDS(expr_gene, file.path(DIR_PROC, "expr_gene.rds"))
readr::write_csv(
  tibble::as_tibble(expr_gene, rownames = "gene"),
  file.path(DIR_PROC, "expr_gene.csv")
)

## ---- 5. Integration QC -----------------------------------------------------
## Gate 1: cross-platform corneal concordance must reach the pre-specified
##         threshold, otherwise the harmonization is not fit for the C1/C3
##         cross-platform contrasts.

r_cornea <- stats::cor(
  rowMeans(expr_gene[, cornea_570,  drop = FALSE]),
  rowMeans(expr_gene[, cornea_6244, drop = FALSE]),
  method = "pearson", use = "pairwise.complete.obs"
)
message(sprintf("cross-platform corneal concordance: r = %.3f", r_cornea))
assert_that(
  r_cornea >= CONCORDANCE_MIN,
  sprintf("Corneal concordance r = %.3f is below the pre-specified %.2f threshold.",
          r_cornea, CONCORDANCE_MIN)
)

## Gate 2: PCA should separate samples by tissue, not by batch.
pca  <- stats::prcomp(t(expr_gene), scale. = TRUE)
vexp <- pca$sdev^2 / sum(pca$sdev^2)

pca_out <- tibble::tibble(
  gsm      = rownames(pca$x),
  group    = meta$group[match(rownames(pca$x), meta$gsm)],
  series   = meta$series[match(rownames(pca$x), meta$gsm)],
  PC1      = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3]
)
attr(pca_out, "var_explained") <- vexp[1:3]
write_result(pca_out, file.path(DIR_PROC, "pca_coordinates.csv"))

write_result(
  tibble::tibble(
    metric = c("n_transcripts", "n_genes", "cornea_concordance_r",
               "PC1_var_explained", "PC2_var_explained", "PC3_var_explained"),
    value  = c(nrow(expr_tx), nrow(expr_gene), round(r_cornea, 4),
               round(vexp[1], 4), round(vexp[2], 4), round(vexp[3], 4))
  ),
  file.path(DIR_DE, "integration_qc.csv")
)

message("02_integration.R complete.")
