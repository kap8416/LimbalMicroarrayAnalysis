## =============================================================================
## 07_core_signature_null.R
## Permutation null for the core limbal signature.
##
## WHY THIS EXISTS
## The three contrasts share the same four limbal samples, so they are not
## independent tests. A hypergeometric or independence-based enrichment for the
## three-way intersection is therefore not interpretable: a gene that happens to
## be high in those particular biopsies tends to appear upregulated in all three
## contrasts by construction.
##
## This script builds the correct null by permuting group labels and recomputing
## the entire three-contrast pipeline, so the null intersection inherits the same
## shared-sample dependency as the observed one.
##
## Labels are permuted WITHIN PLATFORM. Permuting across platforms would break the
## confounding structure and produce an optimistically small null, because the
## permuted contrasts would no longer carry the platform component that inflates
## the observed ones.
##
## Output:
##   results/de/core_signature_null.csv    null intersection sizes
##   results/de/core_signature_null.txt    observed value, null summary, p value
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(limma)
})

N_PERM <- 1000L

expr <- readRDS(file.path(DIR_PROC, "expr_gene.rds"))
meta <- readr::read_csv(file.path(DIR_RAW, "sample_metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(group = factor(group, levels = GROUP_LEVELS))
meta <- meta[match(colnames(expr), meta$gsm), ]

## ---- Core-signature size for a given grouping ------------------------------

#' Number of genes upregulated in the reference group in all three contrasts.
#'
#' Applies the same thresholds, the same no-intercept design and the same
#' eBayes(trend = TRUE) moderation as 03_differential_expression.R, so the null
#' and the observed statistic are computed identically.
core_size <- function(groups, named_only = TRUE) {
  ## A permutation can leave a group too small to fit; report NA rather than
  ## silently returning a misleading zero.
  if (any(table(groups)[GROUP_LEVELS] < 2L)) return(NA_integer_)

  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)

  fit <- limma::lmFit(expr, design)
  cm <- limma::makeContrasts(
    C1 = Limbus - Cornea,
    C2 = Limbus - Conjunctiva,
    C3 = Limbus - Cultured,
    levels = design
  )
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE)

  up <- lapply(CONTRASTS$id, function(id) {
    cc <- CONTRASTS[CONTRASTS$id == id, ]
    tt <- limma::topTable(fit2, coef = id, number = Inf, adjust.method = "BH")
    rownames(tt)[tt$adj.P.Val <= cc$fdr_cut & tt$logFC >= cc$lfc_cut]
  })

  shared <- Reduce(intersect, up)
  if (named_only) shared <- shared[!grepl("^(XM_|XR_|NR_)", shared)]
  length(shared)
}

## ---- Observed ---------------------------------------------------------------

observed <- core_size(meta$group)
message("observed core signature (named genes): ", observed)

## ---- Null distribution -----------------------------------------------------
## Within-platform permutation: samples are shuffled among the groups present on
## their own platform, preserving both group sizes and the platform-group
## confounding.

set.seed(SEED)
plat <- meta$platform

null_sizes <- vapply(seq_len(N_PERM), function(i) {
  g <- meta$group
  for (p in unique(plat)) {
    idx <- which(plat == p)
    g[idx] <- sample(g[idx])
  }
  if (i %% 50L == 0L) message("  permutation ", i, "/", N_PERM)
  core_size(g)
}, integer(1))

null_sizes <- null_sizes[!is.na(null_sizes)]
n_valid <- length(null_sizes)

## Add-one estimator: a permutation p value can never be exactly zero.
p_emp <- (sum(null_sizes >= observed) + 1) / (n_valid + 1)

write_result(
  tibble::tibble(permutation = seq_along(null_sizes), core_size = null_sizes),
  file.path(DIR_DE, "core_signature_null.csv")
)

summary_txt <- c(
  "Permutation null for the core limbal signature",
  "==============================================",
  sprintf("Permutations requested / valid : %d / %d", N_PERM, n_valid),
  sprintf("Permutation scheme             : group labels shuffled within platform"),
  sprintf("Statistic                      : genes upregulated in all three contrasts"),
  sprintf("Thresholds                     : C1 %.1f, C2 %.1f, C3 %.1f; FDR <= 0.05",
          CONTRASTS$lfc_cut[1], CONTRASTS$lfc_cut[2], CONTRASTS$lfc_cut[3]),
  "",
  sprintf("Observed                       : %d", observed),
  sprintf("Null mean (SD)                 : %.2f (%.2f)",
          mean(null_sizes), stats::sd(null_sizes)),
  sprintf("Null median [IQR]              : %.0f [%.0f-%.0f]",
          stats::median(null_sizes),
          stats::quantile(null_sizes, 0.25), stats::quantile(null_sizes, 0.75)),
  sprintf("Null 95th percentile           : %.0f", stats::quantile(null_sizes, 0.95)),
  sprintf("Null max                       : %d", max(null_sizes)),
  sprintf("Permutations >= observed       : %d", sum(null_sizes >= observed)),
  sprintf("Empirical p                    : %.4g%s",
          p_emp, if (sum(null_sizes >= observed) == 0) " (upper bound)" else ""),
  "",
  "Interpretation note: this null preserves the fact that all three contrasts",
  "share the same limbal samples. It therefore tests whether the three-way",
  "intersection is larger than expected given that dependency, which an",
  "independence-based (hypergeometric) calculation does not."
)

writeLines(summary_txt, file.path(DIR_DE, "core_signature_null.txt"))
cat(summary_txt, sep = "\n")

message("\n07_core_signature_null.R complete.")
