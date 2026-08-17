## =============================================================================
## 06_mra_rewiring.R
## Co-expression network inference, master regulator analysis (MRA), and
## condition-specific regulatory rewiring.
##
## MRA scores REGULON ACTIVITY, not transcript abundance. A regulator can
## therefore be significantly repressed by MRA while its own transcript is
## constitutively expressed and non-differential (PAX6 is the case in point).
## This distinction must be preserved in the manuscript wording.
##
## Output:
##   results/mra/regulon_summary.csv        regulon sizes per TF
##   results/mra/MRA_<contrast>.csv         NES, p, FDR per regulator
##   results/mra/MRA_all.csv                concatenated
##   results/mra/rewiring_edges.csv         rewired edges (|delta rho| >= 0.6)
##   results/mra/rewiring_tf_summary.csv    rewired-edge count per TF
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(corto)
})

expr <- readRDS(file.path(DIR_PROC, "expr_gene.rds"))
meta <- readr::read_csv(file.path(DIR_RAW, "sample_metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(group = factor(group, levels = GROUP_LEVELS))
meta <- meta[match(colnames(expr), meta$gsm), ]

## ---- 0. Sample-set gate ----------------------------------------------------
## The original analysis ran this layer on a 17-sample matrix that retained the
## conjunctival outlier GSM724094, while the differential expression layer used
## the curated 16 (docs/AUDIT.md 3.5). Conjunctiva was therefore n = 4 here and
## n = 3 there, and every published NES and rewired-edge count carries that
## inconsistency. These assertions make the divergence impossible to repeat.

assert_that(
  ncol(expr) == 16L,
  sprintf(paste0("Expected the curated 16-sample matrix; got %d columns. Every ",
                 "analytical layer must use the same sample set."), ncol(expr))
)
assert_that(
  !any(EXCLUDED_SAMPLES$gsm %in% colnames(expr)),
  paste0("Excluded samples present in the expression matrix: ",
         paste(intersect(EXCLUDED_SAMPLES$gsm, colnames(expr)), collapse = ", "),
         ". These were removed during curation and must not reappear here.")
)
assert_that(!any(is.na(meta$group)), "Sample metadata does not cover all matrix columns.")

obs_n <- table(meta$group)
assert_that(
  all(obs_n[names(EXPECTED_N)] == EXPECTED_N),
  paste0("Group sizes differ from the curated design. Observed: ",
         paste(names(obs_n), obs_n, sep = "=", collapse = ", "))
)
message("sample-set gate passed: ",
        paste(names(obs_n), obs_n, sep = "=", collapse = ", "))

## ---- 1. Input matrix -------------------------------------------------------
## Restricting to the most variable genes concentrates the mutual-information
## signal that corto uses to build regulons and removes invariant features that
## contribute noise rather than structure.

vars <- apply(expr, 1, stats::var, na.rm = TRUE)
top  <- names(sort(vars, decreasing = TRUE))[seq_len(min(CORTO_NTOP, length(vars)))]
mat  <- expr[top, , drop = FALSE]
message("corto input: ", nrow(mat), " genes x ", ncol(mat), " samples")

## ---- 2. Transcription factor centroids -------------------------------------
## Centroids are the annotated human TFs present in the input matrix. Replace
## the fallback list with a curated resource (e.g. the Lambert et al. human TF
## census) for a publication-grade run; the file is read if present.

tf_file <- file.path(DIR_RAW, "human_tfs.txt")
assert_that(
  file.exists(tf_file),
  paste0("Transcription-factor census not found at ", tf_file, ".\n",
         "Run R/00b_build_tf_census.R first. Master regulator analysis must not be ",
         "run on an ad hoc centroid list: the candidate set determines what the ",
         "analysis can find, and a hand-picked one makes the result circular ",
         "(docs/AUDIT.md 3.9).")
)
tfs <- unique(trimws(readLines(tf_file, warn = FALSE)))
tfs <- tfs[nzchar(tfs)]
message("TF census: ", length(tfs), " factors from ", basename(tf_file))
assert_that(
  length(tfs) >= 800L,
  sprintf(paste0("The census holds only %d factors; a human TF census holds roughly ",
                 "1,600. Re-run R/00b_build_tf_census.R."), length(tfs))
)

centroids <- intersect(tfs, rownames(mat))
message(sprintf("centroids present in the top-%d variable genes: %d",
                CORTO_NTOP, length(centroids)))
## corto needs enough centroids for the regulon inference to be meaningful. If the
## variable-gene filter leaves too few, widen CORTO_NTOP rather than proceeding.
assert_that(
  length(centroids) >= 150L,
  sprintf(paste0("Only %d TF centroids survive the top-%d variable-gene filter. ",
                 "Increase CORTO_NTOP in R/00_config.R."),
          length(centroids), CORTO_NTOP)
)

## ---- 3. Regulon inference --------------------------------------------------

regulon_rds <- file.path(DIR_PROC, "corto_regulon.rds")
if (file.exists(regulon_rds)) {
  regulon <- readRDS(regulon_rds)
  message("cached regulon: ", length(regulon), " regulators")
} else {
  set.seed(SEED)
  regulon <- corto::corto(
    inmat     = mat,
    centroids = centroids,
    nbootstraps = CORTO_NBOOTSTRAP,
    p         = CORTO_PVALUE,
    nthreads  = max(1L, parallel::detectCores() - 1L),
    verbose   = TRUE
  )
  saveRDS(regulon, regulon_rds)
}

regulon_summary <- tibble::tibble(
  tf   = names(regulon),
  size = vapply(regulon, function(r) length(r$tfmode), integer(1))
) |>
  dplyr::arrange(dplyr::desc(size))
write_result(regulon_summary, file.path(DIR_MRA, "regulon_summary.csv"))

## ---- 4. Master regulator analysis per contrast -----------------------------

mra_all <- list()

for (i in seq_len(nrow(CONTRASTS))) {
  cc <- CONTRASTS[i, ]
  s1 <- meta$gsm[meta$group == cc$reference]
  s2 <- meta$gsm[meta$group == cc$comparator]

  message(sprintf("MRA %s: %s (n=%d) vs %s (n=%d)",
                  cc$id, cc$reference, length(s1), cc$comparator, length(s2)))

  res <- corto::mra(
    expmat1 = mat[, s1, drop = FALSE],
    expmat2 = mat[, s2, drop = FALSE],
    regulon = regulon,
    minsize = MRA_MINSIZE,
    nperm   = 1000,
    verbose = FALSE
  )

  tab <- tibble::tibble(
    contrast     = cc$id,
    label        = cc$label,
    tf           = names(res$nes),
    NES          = as.numeric(res$nes),
    pvalue       = as.numeric(res$pvalue[names(res$nes)]),
    regulon_size = regulon_summary$size[match(names(res$nes), regulon_summary$tf)]
  ) |>
    dplyr::mutate(
      padj      = stats::p.adjust(pvalue, method = "BH"),
      activity  = dplyr::if_else(NES > 0, "Activated in limbus", "Repressed in limbus"),
      ## Explicit reminder that NES reflects regulon activity, not TF abundance.
      interpretation = "NES = regulon activity; not transcript abundance"
    ) |>
    dplyr::arrange(dplyr::desc(abs(NES)))

  write_result(tab, file.path(DIR_MRA, sprintf("MRA_%s.csv", cc$id)))
  mra_all[[cc$id]] <- tab

  top <- tab |> dplyr::slice_head(n = 5)
  message("  top |NES|: ",
          paste(sprintf("%s (%.2f)", top$tf, top$NES), collapse = ", "))
}

write_result(dplyr::bind_rows(mra_all), file.path(DIR_MRA, "MRA_all.csv"))

## ---- 5. Condition-specific co-expression rewiring --------------------------
## For each TF-target pair present in the regulons, Spearman correlation is
## recomputed within each group. Edges whose correlation changes by at least
## REWIRE_DELTA_RHO between two conditions are classified as rewired.
##
## Rewiring is deliberately independent of differential expression: a regulator
## can retain constant expression while its target coupling is reorganized, and
## that reorganization is the quantity of interest here.

group_mats <- lapply(setNames(GROUP_LEVELS, GROUP_LEVELS), function(g) {
  mat[, meta$gsm[meta$group == g], drop = FALSE]
})

edge_list <- purrr::imap_dfr(regulon, function(r, tf) {
  tibble::tibble(tf = tf, target = names(r$tfmode))
}) |>
  dplyr::filter(target %in% rownames(mat), tf %in% rownames(mat))

message("evaluating rewiring across ", nrow(edge_list), " regulon edges")

rho_by_group <- vapply(names(group_mats), function(g) {
  m <- group_mats[[g]]
  suppressWarnings(
    mapply(function(a, b) stats::cor(m[a, ], m[b, ], method = "spearman"),
           edge_list$tf, edge_list$target)
  )
}, numeric(nrow(edge_list)))

rewiring <- dplyr::bind_cols(edge_list, tibble::as_tibble(rho_by_group)) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    delta_max = max(dplyr::c_across(dplyr::all_of(GROUP_LEVELS)), na.rm = TRUE) -
                min(dplyr::c_across(dplyr::all_of(GROUP_LEVELS)), na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(rewired = delta_max >= REWIRE_DELTA_RHO)

write_result(dplyr::filter(rewiring, rewired), file.path(DIR_MRA, "rewiring_edges.csv"))

rewiring_tf <- rewiring |>
  dplyr::filter(rewired) |>
  dplyr::count(tf, name = "n_rewired_edges") |>
  dplyr::left_join(regulon_summary, by = "tf") |>
  dplyr::mutate(frac_rewired = round(n_rewired_edges / size, 3)) |>
  dplyr::arrange(dplyr::desc(n_rewired_edges))
write_result(rewiring_tf, file.path(DIR_MRA, "rewiring_tf_summary.csv"))

message("most rewired regulators: ",
        paste(sprintf("%s (%d)", rewiring_tf$tf[1:9], rewiring_tf$n_rewired_edges[1:9]),
              collapse = ", "))

message("06_mra_rewiring.R complete.")
