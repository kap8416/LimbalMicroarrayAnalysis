## =============================================================================
## 02b_use_existing_matrix.R
## Alternative to 02_integration.R: start from the harmonized matrix produced by
## the original analysis instead of regenerating it.
##
## WHY THIS EXISTS
## The script that produced combined_expression_matrix_corrected_v2.csv was never
## recovered (docs/AUDIT.md 3.7). Re-running 02_integration.R would regenerate the
## matrix from the Methods description, but the result would not be byte-identical
## to the one the published tables came from, so any change in downstream numbers
## would be impossible to attribute.
##
## Using this script instead isolates the problem: everything downstream is
## re-run with the corrected code, on exactly the input the original analysis used.
## Differences in the results can then be attributed to the corrections rather than
## to a re-derived matrix.
##
## Run EITHER 02_integration.R OR this script, not both.
##
## Input:
##   data/raw/combined_expression_matrix_corrected_v2.csv   (45,378 x 20, RefSeq IDs)
##
## Output:
##   data/processed/expr_transcript.rds   curated 16 samples, transcript level
##   data/processed/expr_gene.rds         gene level, highest-MAD isoform per symbol
##   data/processed/expr_gene.csv
##   data/processed/pca_coordinates.csv
##   results/de/integration_qc.csv
##   results/de/duplicate_profile_report.csv
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

MATRIX_FILE <- file.path(DIR_RAW, "combined_expression_matrix_corrected_v2.csv")
assert_that(
  file.exists(MATRIX_FILE),
  paste0("Not found: ", MATRIX_FILE,
         "\nCopy the original harmonized matrix there, or run R/02_integration.R ",
         "instead to regenerate it from GEO.")
)

meta <- readr::read_csv(file.path(DIR_RAW, "sample_metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(group = factor(group, levels = GROUP_LEVELS))

## ---- 1. Load and curate ----------------------------------------------------

raw <- readr::read_csv(MATRIX_FILE, show_col_types = FALSE)
ids <- raw[[1]]
expr_all <- as.matrix(raw[, -1])
rownames(expr_all) <- ids
message(sprintf("loaded %d transcripts x %d samples", nrow(expr_all), ncol(expr_all)))

## The file carries all 20 samples; restrict to the curated 16. This is the step
## the original corto script omitted, which is how the outlier reached the
## regulatory layer (docs/AUDIT.md 3.5).
assert_that(
  all(meta$gsm %in% colnames(expr_all)),
  paste0("Curated samples missing from the matrix: ",
         paste(setdiff(meta$gsm, colnames(expr_all)), collapse = ", "))
)
dropped <- setdiff(colnames(expr_all), meta$gsm)
message("dropping ", length(dropped), " non-analytical samples: ",
        paste(dropped, collapse = ", "))

expr_tx <- expr_all[, meta$gsm, drop = FALSE]
assert_that(ncol(expr_tx) == 16L, "Expected 16 curated samples.")
saveRDS(expr_tx, file.path(DIR_PROC, "expr_transcript.rds"))

## ---- 2. Duplicate-profile audit -------------------------------------------
## 82.7% of the rows in this matrix are exact duplicates of another row: many
## RefSeq accessions are transcript variants measured by the same probe, so their
## values were copied rather than measured independently.
##
## This matters downstream. Distinct gene symbols that inherit identical values
## are not independent observations: they enter DEG lists together and, if they
## share annotation terms, they inflate over-representation significance. The
## report written here quantifies the effect so it can be stated rather than
## discovered by a reviewer.

prof <- apply(round(expr_tx, 10), 1, paste, collapse = "|")
dup_tx <- sum(duplicated(prof) | duplicated(prof, fromLast = TRUE))
message(sprintf("transcript level: %d/%d rows share an identical profile (%.1f%%); %d unique profiles",
                dup_tx, nrow(expr_tx), 100 * dup_tx / nrow(expr_tx),
                length(unique(prof))))

## ---- 3. Transcript -> gene collapse by maximum MAD -------------------------
## RefSeq accession -> symbol via org.Hs.eg.db, so this step does not depend on a
## platform annotation file. Accessions without a symbol keep their accession as
## identifier, which is why predicted XM_/XR_ features persist.

acc <- sub("\\..*$", "", rownames(expr_tx))   # strip version suffixes if present
map <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = acc, column = "SYMBOL", keytype = "REFSEQ", multiVals = "first"
))
sym <- unname(map)
sym[is.na(sym)] <- rownames(expr_tx)[is.na(sym)]
message(sprintf("mapped %d/%d accessions to a gene symbol",
                sum(!is.na(map)), length(acc)))

mads <- row_mad(expr_tx)
## Ties (identical isoform profiles) resolve to the first accession in input
## order, which is deterministic and harmless because the values are identical.
keep <- tapply(seq_len(nrow(expr_tx)), sym, function(i) i[which.max(mads[i])])
expr_gene <- expr_tx[unlist(keep), , drop = FALSE]
rownames(expr_gene) <- names(keep)

message(sprintf("gene level: %d features", nrow(expr_gene)))
if (nrow(expr_gene) != N_GENES_EXPECTED) {
  warning(sprintf(paste0(
    "Gene count is %d; the published analysis reports %d. The difference comes ",
    "from the symbol mapping resource (org.Hs.eg.db here versus the GPL570 ",
    "annotation table originally). Record whichever is used."),
    nrow(expr_gene), N_GENES_EXPECTED))
}

## How much duplication survives the collapse?
prof_g <- apply(round(expr_gene, 10), 1, paste, collapse = "|")
dup_g  <- sum(duplicated(prof_g) | duplicated(prof_g, fromLast = TRUE))
message(sprintf("gene level: %d/%d features share an identical profile (%.1f%%); %d unique profiles",
                dup_g, nrow(expr_gene), 100 * dup_g / nrow(expr_gene),
                length(unique(prof_g))))

write_result(
  tibble::tibble(
    level            = c("transcript", "gene"),
    n_features       = c(nrow(expr_tx), nrow(expr_gene)),
    n_unique_profiles = c(length(unique(prof)), length(unique(prof_g))),
    n_duplicated     = c(dup_tx, dup_g),
    pct_duplicated   = round(100 * c(dup_tx / nrow(expr_tx), dup_g / nrow(expr_gene)), 1)
  ),
  file.path(DIR_DE, "duplicate_profile_report.csv")
)

saveRDS(expr_gene, file.path(DIR_PROC, "expr_gene.rds"))
readr::write_csv(tibble::as_tibble(expr_gene, rownames = "gene"),
                 file.path(DIR_PROC, "expr_gene.csv"))

## ---- 4. Integration QC on the gene-level matrix ---------------------------

cornea_570  <- meta$gsm[meta$platform == "GPL570"  & meta$group == "Cornea"]
cornea_6244 <- meta$gsm[meta$platform == "GPL6244" & meta$group == "Cornea"]

r_cornea <- stats::cor(
  rowMeans(expr_gene[, cornea_570,  drop = FALSE]),
  rowMeans(expr_gene[, cornea_6244, drop = FALSE]),
  method = "pearson", use = "pairwise.complete.obs"
)
message(sprintf("cross-platform corneal concordance (gene level): r = %.4f", r_cornea))
assert_that(
  r_cornea >= CONCORDANCE_MIN,
  sprintf("Corneal concordance r = %.4f is below the pre-specified %.2f.",
          r_cornea, CONCORDANCE_MIN)
)

pca  <- stats::prcomp(t(expr_gene), scale. = TRUE)
vexp <- pca$sdev^2 / sum(pca$sdev^2)

pca_out <- tibble::tibble(
  gsm      = rownames(pca$x),
  group    = meta$group[match(rownames(pca$x), meta$gsm)],
  platform = meta$platform[match(rownames(pca$x), meta$gsm)],
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3]
)
write_result(pca_out, file.path(DIR_PROC, "pca_coordinates.csv"))

## Does any leading PC separate samples by platform rather than by tissue? A large
## value here would mean residual batch structure survived correction.
plat_sep <- vapply(1:3, function(k) {
  a <- pca$x[pca_out$platform == "GPL570",  k]
  b <- pca$x[pca_out$platform == "GPL6244", k]
  abs(mean(a) - mean(b)) / stats::sd(pca$x[, k])
}, numeric(1))
message(sprintf("platform separation: PC1 %.2f, PC2 %.2f, PC3 %.2f (SD units)",
                plat_sep[1], plat_sep[2], plat_sep[3]))

## Within-group sample correlations, the analysis that justified the exclusion.
cors <- stats::cor(expr_gene, method = "pearson")
within <- purrr::map_dfr(levels(meta$group), function(g) {
  s <- meta$gsm[meta$group == g]
  if (length(s) < 2) return(NULL)
  purrr::map_dfr(s, function(x) {
    o <- setdiff(s, x)
    tibble::tibble(gsm = x, group = g,
                   mean_r_within = mean(cors[x, o]),
                   min_r_within  = min(cors[x, o]))
  })
})
write_result(within, file.path(DIR_DE, "within_group_correlation.csv"))

write_result(
  tibble::tibble(
    metric = c("n_transcripts", "n_genes", "n_unique_gene_profiles",
               "cornea_concordance_r",
               "PC1_var_explained", "PC2_var_explained", "PC3_var_explained",
               "platform_sep_PC1", "platform_sep_PC2", "platform_sep_PC3",
               "matrix_source"),
    value  = c(nrow(expr_tx), nrow(expr_gene), length(unique(prof_g)),
               round(r_cornea, 4),
               round(vexp[1], 4), round(vexp[2], 4), round(vexp[3], 4),
               round(plat_sep[1], 3), round(plat_sep[2], 3), round(plat_sep[3], 3),
               basename(MATRIX_FILE))
  ),
  file.path(DIR_DE, "integration_qc.csv")
)

message("\n02b_use_existing_matrix.R complete.")
message("Downstream steps (03 onwards) can now run unchanged.")
