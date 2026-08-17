## =============================================================================
## 04_functional_enrichment.R
## GO Biological Process and KEGG over-representation analysis (ORA) with
## clusterProfiler, run separately for up- and downregulated gene sets in each
## contrast.
##
## The universe is restricted to genes present in the batch-corrected matrix.
## Using the full genome as background would inflate significance, because the
## measurable gene space is bounded by the intersection of the two platforms.
##
## Output:
##   results/enrichment/GO_BP_<contrast>_<direction>.csv
##   results/enrichment/KEGG_<contrast>_<direction>.csv
##   results/enrichment/enrichment_all.csv    concatenated, tidy
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

expr       <- readRDS(file.path(DIR_PROC, "expr_gene.rds"))
membership <- readr::read_csv(file.path(DIR_DE, "deg_membership.csv"), show_col_types = FALSE)

## ---- 1. Symbol -> Entrez mapping ------------------------------------------
## Predicted accessions carry no Entrez ID and drop out here, which is expected.

to_entrez <- function(symbols) {
  suppressWarnings(
    clusterProfiler::bitr(unique(symbols), fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )$ENTREZID |> unique()
}

universe_entrez <- to_entrez(rownames(expr))
message("universe: ", length(universe_entrez), " Entrez IDs")

## ---- 2. ORA per contrast and direction ------------------------------------

run_ora <- function(gene_entrez, label) {
  if (length(gene_entrez) < 5L) {
    message("  ", label, ": fewer than 5 mapped genes, skipped")
    return(list(go = NULL, kegg = NULL))
  }

  go <- clusterProfiler::enrichGO(
    gene          = gene_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = GO_ONTOLOGY,
    keyType       = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff  = ENR_FDR_CUT,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )

  kegg <- clusterProfiler::enrichKEGG(
    gene          = gene_entrez,
    universe      = universe_entrez,
    organism      = "hsa",
    keyType       = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff  = ENR_FDR_CUT,
    qvalueCutoff  = 0.20
  )
  if (!is.null(kegg)) kegg <- clusterProfiler::setReadable(kegg, org.Hs.eg.db, "ENTREZID")

  list(go = go, kegg = kegg)
}

tidy_enrich <- function(obj, contrast, direction, source) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0L) return(NULL)
  as.data.frame(obj) |>
    tibble::as_tibble() |>
    dplyr::transmute(
      contrast, direction, source,
      term_id = ID, term = Description,
      gene_ratio = GeneRatio, bg_ratio = BgRatio,
      pvalue, padj = p.adjust, qvalue,
      count = Count, genes = geneID
    )
}

all_enrich <- list()

for (cid in CONTRASTS$id) {
  status_col <- paste0(cid, "_status")

  for (dir_lab in c("Up", "Down")) {
    genes  <- membership$Symbol[membership[[status_col]] == dir_lab]
    genes  <- genes[!is.na(genes)]
    entrez <- to_entrez(genes)
    label  <- paste0(cid, " ", dir_lab)
    message(sprintf("ORA %s: %d genes -> %d Entrez", label, length(genes), length(entrez)))

    res <- run_ora(entrez, label)

    go_t   <- tidy_enrich(res$go,   cid, dir_lab, "GO_BP")
    kegg_t <- tidy_enrich(res$kegg, cid, dir_lab, "KEGG")

    if (!is.null(go_t))   write_result(go_t,   file.path(DIR_ENR, sprintf("GO_BP_%s_%s.csv", cid, dir_lab)))
    if (!is.null(kegg_t)) write_result(kegg_t, file.path(DIR_ENR, sprintf("KEGG_%s_%s.csv",  cid, dir_lab)))

    all_enrich <- c(all_enrich, list(go_t, kegg_t))
  }
}

enrich_all <- dplyr::bind_rows(all_enrich)

## ---- 3. Interpretation flags ----------------------------------------------
## Two GO terms recur in these contrasts with names that invite over-reading.
## They are flagged here so the annotation, rather than the term label, drives
## interpretation downstream:
##   - "ossification"          : driven by shared ECM and BMP-pathway genes, not
##                               by osteogenic differentiation of the limbus.
##   - "epidermis development" : reflects keratinization programs common to
##                               stratified epithelia, not epidermal specification.
FLAGGED_TERMS <- c("ossification", "epidermis development", "keratinocyte differentiation")

enrich_all <- enrich_all |>
  dplyr::mutate(
    needs_caveat = tolower(term) %in% FLAGGED_TERMS,
    caveat = dplyr::case_when(
      tolower(term) == "ossification" ~
        "Annotation driven by shared ECM/BMP-pathway genes, not osteogenic differentiation.",
      tolower(term) %in% c("epidermis development", "keratinocyte differentiation") ~
        "GO term name reflects keratinization programs shared by stratified epithelia, not epidermal specification.",
      TRUE ~ NA_character_
    )
  )

write_result(enrich_all, file.path(DIR_ENR, "enrichment_all.csv"))

if (any(enrich_all$needs_caveat)) {
  message("\nTerms requiring an interpretation caveat in the manuscript text:")
  enrich_all |>
    dplyr::filter(needs_caveat) |>
    dplyr::select(contrast, direction, term, padj, count) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

message("04_functional_enrichment.R complete.")
