## =============================================================================
## 03_differential_expression.R
## limma differential expression for the three limbus-anchored contrasts, plus
## derivation of the core limbal signature shared across all three.
##
## eBayes(trend = TRUE) applies an intensity-dependent prior variance, which is
## the appropriate empirical-Bayes moderation for normalized log2 microarray
## intensities. voom is NOT used here: it models the mean-variance relationship
## of counts and is not applicable to microarray data.
##
## Output:
##   results/de/DEA_<contrast>.csv        full topTable per contrast
##   results/de/DEA_<contrast>_UP.csv     significant, higher in limbus
##   results/de/DEA_<contrast>_DOWN.csv   significant, lower in limbus
##   results/de/DEA_counts_summary.csv    up / down / ns per contrast
##   results/de/core_signature.csv        genes upregulated in limbus in C1-C3
##   results/de/deg_membership.csv        per-gene status and logFC in all three
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(limma)
})

expr <- readRDS(file.path(DIR_PROC, "expr_gene.rds"))
meta <- readr::read_csv(file.path(DIR_RAW, "sample_metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(group = factor(group, levels = GROUP_LEVELS))
meta <- meta[match(colnames(expr), meta$gsm), ]
assert_that(!any(is.na(meta$group)), "Sample metadata does not cover all matrix columns.")

## ---- 1. Fit the linear model ----------------------------------------------
## A no-intercept design gives one coefficient per group, so every contrast is
## expressed directly as a difference of group means.

design <- stats::model.matrix(~ 0 + meta$group)
colnames(design) <- levels(meta$group)

fit <- limma::lmFit(expr, design)

cm <- limma::makeContrasts(
  C1 = Limbus - Cornea,
  C2 = Limbus - Conjunctiva,
  C3 = Limbus - Cultured,
  levels = design
)

fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE)

## ---- 2. Extract per-contrast results --------------------------------------

de_tables <- list()

for (i in seq_len(nrow(CONTRASTS))) {
  cc  <- CONTRASTS[i, ]
  tt  <- limma::topTable(fit2, coef = cc$id, number = Inf, adjust.method = "BH")
  tt  <- tibble::as_tibble(tt, rownames = "Symbol")

  tt <- tt |>
    dplyr::mutate(
      status = dplyr::case_when(
        adj.P.Val <= cc$fdr_cut &  logFC >=  cc$lfc_cut ~ "Up",
        adj.P.Val <= cc$fdr_cut &  logFC <= -cc$lfc_cut ~ "Down",
        TRUE                                            ~ "NS"
      ),
      contrast = cc$id
    )

  de_tables[[cc$id]] <- tt

  stem <- sprintf("DEA_%s_%s_vs_%s", cc$id, cc$reference, cc$comparator)
  write_result(tt,                                    file.path(DIR_DE, paste0(stem, ".csv")))
  write_result(dplyr::filter(tt, status == "Up"),     file.path(DIR_DE, paste0(stem, "_UP.csv")))
  write_result(dplyr::filter(tt, status == "Down"),   file.path(DIR_DE, paste0(stem, "_DOWN.csv")))

  message(sprintf("  %s (%s): %d up, %d down  [|log2FC| >= %.1f, FDR <= %.2f]%s",
                  cc$id, cc$label,
                  sum(tt$status == "Up"), sum(tt$status == "Down"),
                  cc$lfc_cut, cc$fdr_cut,
                  if (cc$exploratory) "  (exploratory)" else ""))
}

## ---- 3. Counts summary -----------------------------------------------------

counts_summary <- purrr::imap_dfr(de_tables, function(tt, id) {
  cc <- CONTRASTS[CONTRASTS$id == id, ]
  tibble::tibble(
    contrast = id,
    label    = cc$label,
    lfc_cut  = cc$lfc_cut,
    fdr_cut  = cc$fdr_cut,
    up       = sum(tt$status == "Up"),
    down     = sum(tt$status == "Down"),
    ns       = sum(tt$status == "NS"),
    total    = nrow(tt)
  )
})
write_result(counts_summary, file.path(DIR_DE, "DEA_counts_summary.csv"))

## ---- 4. Cross-contrast membership -----------------------------------------
## One row per gene that is differentially expressed in at least one contrast,
## with status and effect size in all three. This is the table consumed by the
## overlap figure and by the core-signature derivation.

membership <- purrr::reduce(
  purrr::imap(de_tables, function(tt, id) {
    tt |>
      dplyr::select(Symbol, logFC, adj.P.Val, status) |>
      dplyr::rename_with(~ paste0(id, "_", .x), -Symbol)
  }),
  dplyr::full_join, by = "Symbol"
)

membership <- membership |>
  dplyr::mutate(
    n_up = (C1_status == "Up") + (C2_status == "Up") + (C3_status == "Up"),
    module = dplyr::case_when(
      C1_status == "Up"   & C2_status == "Up"   & C3_status == "Up"   ~ "CORE - all 3 contrasts",
      C1_status == "Up"   & C2_status == "Up"                        ~ "SHARED C1 & C2",
      C1_status == "Up"   & C3_status == "Up"                        ~ "SHARED C1 & C3",
      C2_status == "Up"   & C3_status == "Up"                        ~ "SHARED C2 & C3",
      C1_status == "Up"                                             ~ "C1 SPECIFIC",
      C2_status == "Up"                                             ~ "C2 SPECIFIC",
      C3_status == "Up"                                             ~ "C3 SPECIFIC",
      TRUE                                                          ~ "DOWN IN LIMBUS"
    )
  ) |>
  dplyr::filter(C1_status != "NS" | C2_status != "NS" | C3_status != "NS")

write_result(membership, file.path(DIR_DE, "deg_membership.csv"))

## ---- 5. Core limbal signature ---------------------------------------------
## Requiring consistent upregulation across three biologically independent
## transitions filters for programs intrinsic to the limbal state rather than
## artefacts of any single comparison.
##
## Predicted-only accessions (XM_/XR_) are reported separately: they satisfy the
## statistical criterion but lack an approved gene symbol, and the manuscript
## reports the named-gene signature.

core_all <- membership |>
  dplyr::filter(module == "CORE - all 3 contrasts") |>
  dplyr::arrange(dplyr::desc(C1_logFC))

core_all <- core_all |>
  dplyr::mutate(is_predicted = grepl("^(XM_|XR_|NR_)", Symbol))

core_named <- dplyr::filter(core_all, !is_predicted)

write_result(core_all,   file.path(DIR_DE, "core_signature_all_features.csv"))
write_result(core_named, file.path(DIR_DE, "core_signature.csv"))

message(sprintf("core signature: %d named genes (+%d predicted transcripts)",
                nrow(core_named), sum(core_all$is_predicted)))

saveRDS(fit2, file.path(DIR_PROC, "limma_fit.rds"))
message("03_differential_expression.R complete.")
