## =============================================================================
## 05_ppi_networks.R
## Contrast-specific protein-protein interaction networks from STRING v12.0,
## with igraph topology, composite hub scoring, and Louvain community detection.
##
## -----------------------------------------------------------------------------
## TWO CORRECTIONS RELATIVE TO THE ORIGINAL ANALYSIS
## (see docs/AUDIT.md sections 3.3 and 3.4)
## -----------------------------------------------------------------------------
##
## 1. COMPLETE INTERACTION RETRIEVAL, NEVER CHUNKED.
##    The original pipeline split the gene list into disjoint chunks of 1,500 and
##    queried the STRING REST API once per chunk. That endpoint returns only
##    interactions AMONG THE SUBMITTED IDENTIFIERS, so every interaction spanning
##    two chunks was never queried. The C1 and C3 networks lost 12% and 5% of
##    their edges, and all reported topology was computed on the incomplete graphs.
##
##    This script therefore prefers the local STRING flat files, which contain the
##    complete interactome and have no query-size limit. The REST API is used only
##    when the whole gene set fits in a single request; if it does not, the script
##    STOPS with instructions rather than degrading silently.
##
## 2. TRUE MAXIMUM CLIQUE CENTRALITY, OR NONE AT ALL.
##    The original pipeline substituted degree^2 for MCC whenever n > 500, which
##    made two of the five "independent" composite terms functions of degree alone.
##    Here MCC is computed exactly. If the clique enumeration exceeds an explicit
##    budget, MCC is dropped from the composite and the omission is recorded in the
##    output — a degree surrogate is never silently substituted.
##
## Output:
##   results/ppi/edges_<contrast>.csv        edge list of the retained component
##   results/ppi/centrality_<contrast>.csv   per-node topology, module, hub score
##   results/ppi/network_summary.csv         nodes, edges, modules, %up, provenance
##   results/ppi/module_summary_<contrast>.csv
## =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(igraph)
})

membership <- readr::read_csv(file.path(DIR_DE, "deg_membership.csv"), show_col_types = FALSE)

## ---- 1. Interaction source -------------------------------------------------

#' Load the complete STRING interactome from the local flat files.
#'
#' Requires two files in data/raw/ (download once, ~90 MB combined):
#'   https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
#'   https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz
#'
#' Returns a data frame of symbol-symbol edges with combined scores on the 0-1
#' scale, filtered to STRING_SCORE_MIN. Cached as an RDS after the first parse.
string_local_edges <- function() {
  cache <- file.path(DIR_PROC, sprintf("string_local_v%s_score%d.rds",
                                       STRING_VERSION, STRING_SCORE_MIN))
  if (file.exists(cache)) {
    message("  cached local STRING interactome: ", basename(cache))
    return(readRDS(cache))
  }

  links_f   <- file.path(DIR_RAW, STRING_LOCAL_LINKS)
  aliases_f <- file.path(DIR_RAW, STRING_LOCAL_ALIASES)
  assert_that(
    file.exists(links_f) && file.exists(aliases_f),
    paste0(
      "Local STRING files not found. Download them once into data/raw/:\n",
      "  curl -LO https://stringdb-downloads.org/download/protein.links.v",
      STRING_VERSION, "/", STRING_LOCAL_LINKS, "\n",
      "  curl -LO https://stringdb-downloads.org/download/protein.aliases.v",
      STRING_VERSION, "/", STRING_LOCAL_ALIASES, "\n",
      "Alternatively set STRING_SOURCE <- \"api\" in R/00_config.R, which is only ",
      "valid if every contrast's gene set fits in a single API request."
    )
  )

  message("  parsing local STRING aliases ...")
  aliases <- readr::read_tsv(
    aliases_f, col_names = c("string_id", "alias", "source"),
    skip = 1, show_col_types = FALSE, progress = FALSE
  )

  ## Sources are matched EXACTLY, in priority order.
  ##
  ## Substring matching is wrong here and fails silently. The v12 aliases file
  ## contains twenty sources beginning "Ensembl_HGNC", of which only
  ## "Ensembl_HGNC_symbol" holds the gene symbol; "Ensembl_HGNC_trans_name" holds
  ## transcript names such as "FN1-201", and others hold CCDS, OMIM, UniProt and
  ## RefSeq identifiers. A pattern like "Ensembl_HGNC" therefore matches all of
  ## them, and deduplicating by file order assigns many proteins a transcript name
  ## instead of a symbol. Nothing then joins to the gene-level DEG set, and the
  ## failure surfaces only as an empty network.
  SYMBOL_SOURCES <- c("Ensembl_HGNC_symbol", "BioMart_HUGO", "Ensembl_EntrezGene",
                      "UniProt_GN_Name", "KEGG_NAME")

  id2sym <- aliases |>
    dplyr::filter(source %in% SYMBOL_SOURCES) |>
    dplyr::mutate(priority = match(source, SYMBOL_SOURCES)) |>
    dplyr::arrange(string_id, priority) |>
    dplyr::distinct(string_id, .keep_all = TRUE)

  assert_that(
    nrow(id2sym) > 15000L,
    sprintf(paste0("Only %d STRING proteins received a gene symbol; expected ~19,000. ",
                   "The aliases file layout may have changed. Sources present: %s"),
            nrow(id2sym), paste(unique(aliases$source)[1:20], collapse = ", "))
  )
  message(sprintf("  %s proteins mapped to a gene symbol", format(nrow(id2sym), big.mark = ",")))

  ## Sanity probe on a gene that must be present and must map to its own symbol.
  fn1 <- id2sym$string_id[id2sym$alias == "FN1"]
  assert_that(length(fn1) >= 1L,
              "FN1 is absent from the STRING symbol map; the mapping is wrong.")

  sym_of <- setNames(id2sym$alias, id2sym$string_id)

  message("  parsing local STRING links (this takes a minute) ...")
  links <- readr::read_delim(
    links_f, delim = " ", show_col_types = FALSE, progress = FALSE
  )
  names(links) <- c("protein1", "protein2", "combined_score")

  edges <- links |>
    dplyr::filter(combined_score >= STRING_SCORE_MIN) |>
    dplyr::transmute(
      from  = sym_of[protein1],
      to    = sym_of[protein2],
      score = combined_score / 1000
    ) |>
    dplyr::filter(!is.na(from), !is.na(to), from != to)

  ## Collapse the reciprocal pairs STRING stores for every interaction.
  edges <- edges |>
    dplyr::mutate(a = pmin(from, to), b = pmax(from, to)) |>
    dplyr::group_by(a, b) |>
    dplyr::summarise(score = max(score), .groups = "drop") |>
    dplyr::rename(from = a, to = b)

  saveRDS(edges, cache)
  message(sprintf("  local interactome: %s interactions at score >= %.2f",
                  format(nrow(edges), big.mark = ","), STRING_SCORE_MIN / 1000))
  edges
}

#' Subset the local interactome to interactions within one gene set.
#'
#' Reports how many of the query genes exist in STRING at all, so a small network
#' can be distinguished from a broken identifier mapping. An empty result with a
#' healthy overlap means the genes are genuinely unconnected; an empty result with
#' no overlap means the identifiers do not match and the run must stop.
string_edges_local <- function(symbols, contrast_id = "") {
  all_edges <- string_local_edges()
  keep <- unique(symbols[!is.na(symbols) & symbols != ""])

  in_string <- sum(keep %in% unique(c(all_edges$from, all_edges$to)))
  message(sprintf("  %d/%d query genes present in the STRING interactome (%.0f%%)",
                  in_string, length(keep), 100 * in_string / max(1L, length(keep))))
  assert_that(
    in_string >= 0.10 * length(keep),
    sprintf(paste0("%s: only %d of %d gene symbols were found in STRING. That is an ",
                   "identifier-mapping failure, not a biological result. Check that the ",
                   "DEG table uses HGNC symbols and that the aliases file parsed correctly."),
            contrast_id, in_string, length(keep))
  )

  dplyr::filter(all_edges, from %in% keep, to %in% keep)
}

#' Single-request STRING REST query.
#'
#' Deliberately NOT chunked. The /network endpoint returns only interactions among
#' the submitted identifiers, so splitting the list discards every cross-chunk
#' interaction. If the set is too large for one request the function stops.
string_edges_api <- function(symbols, contrast_id) {
  cache <- file.path(DIR_RAW, sprintf("string_api_%s_v%s.rds", contrast_id, STRING_VERSION))
  if (file.exists(cache)) {
    message("  cached STRING API response: ", basename(cache))
    return(readRDS(cache))
  }
  assert_that(!OFFLINE,
              sprintf("OFFLINE = TRUE but no cached STRING response for %s", contrast_id))

  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  assert_that(
    length(symbols) <= STRING_API_MAX,
    sprintf(paste0(
      "%s has %d genes, above the %d-identifier single-request limit.\n",
      "Chunking is NOT an option: the STRING /network endpoint returns only ",
      "interactions among the submitted identifiers, so a chunked query silently ",
      "discards every interaction spanning two chunks (this is what produced the ",
      "incomplete C1 and C3 networks; see docs/AUDIT.md 3.3).\n",
      "Set STRING_SOURCE <- \"local\" in R/00_config.R and download the flat files."),
      contrast_id, length(symbols), STRING_API_MAX)
  )

  resp <- httr::RETRY(
    "POST",
    url = paste0(STRING_API, "/json/network"),
    httr::timeout(180),
    body = list(
      identifiers     = paste(symbols, collapse = "\r"),
      species         = STRING_SPECIES,
      required_score  = as.character(STRING_SCORE_MIN),
      network_type    = STRING_NETWORK,
      caller_identity = "LimbalMicroarrayAnalysis"
    ),
    encode = "form",
    times = 3, pause_base = 2
  )
  httr::stop_for_status(resp)

  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  assert_that(nchar(txt) > 3L, sprintf("%s: empty STRING response.", contrast_id))

  edges <- tibble::as_tibble(jsonlite::fromJSON(txt, flatten = TRUE)) |>
    dplyr::transmute(from = preferredName_A, to = preferredName_B, score = score) |>
    dplyr::filter(from != to) |>
    dplyr::mutate(a = pmin(from, to), b = pmax(from, to)) |>
    dplyr::group_by(a, b) |>
    dplyr::summarise(score = max(score), .groups = "drop") |>
    dplyr::rename(from = a, to = b)

  saveRDS(edges, cache)
  edges
}

get_edges <- function(symbols, contrast_id) {
  switch(
    STRING_SOURCE,
    local = string_edges_local(symbols, contrast_id),
    api   = string_edges_api(symbols, contrast_id),
    stop("STRING_SOURCE must be \"local\" or \"api\".", call. = FALSE)
  )
}

## ---- 2. Maximum clique centrality -----------------------------------------

#' Maximum clique centrality (MCC), as defined for cytoHubba.
#'
#'   MCC(v) = sum over maximal cliques C containing v of (|C| - 1)!
#'
#' Factorials are accumulated in log space and exponentiated once, so large
#' cliques do not overflow.
#'
#' Returns NULL if the number of maximal cliques exceeds `budget`. The caller then
#' drops MCC from the composite score rather than substituting a degree surrogate,
#' because degree^2 is a monotone function of degree and would double-weight a
#' term already present in the composite (docs/AUDIT.md 3.4).
mcc_score <- function(g, budget = MCC_CLIQUE_BUDGET) {
  n_cl <- igraph::count_max_cliques(g, min = 2L)
  if (n_cl > budget) {
    message(sprintf("  MCC: %s maximal cliques exceeds the budget of %s; MCC omitted",
                    format(n_cl, big.mark = ","), format(budget, big.mark = ",")))
    return(NULL)
  }
  message(sprintf("  MCC: enumerating %s maximal cliques ...", format(n_cl, big.mark = ",")))

  cl  <- igraph::max_cliques(g, min = 2L)
  out <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)
  for (c in cl) {
    w  <- exp(lgamma(length(c)))          # (|C| - 1)!
    nm <- igraph::V(g)$name[as.integer(c)]
    out[nm] <- out[nm] + w
  }
  ## Nodes in no clique of size >= 2 are isolated within the retained component,
  ## which cannot happen after LCC extraction, but guard anyway.
  out[!is.finite(out)] <- 0
  out
}

minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

## ---- 3. Topology and hub scoring ------------------------------------------

analyse_network <- function(edges, lfc_map, contrast_id) {
  assert_that(nrow(edges) > 0L, sprintf("%s: no interactions retrieved.", contrast_id))

  g <- igraph::graph_from_data_frame(
    dplyr::select(edges, from, to, score), directed = FALSE
  )
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                        edge.attr.comb = list(score = "max"))

  ## Retain the largest connected component: betweenness and closeness are not
  ## comparable across disconnected components.
  comp    <- igraph::components(g)
  n_all   <- igraph::vcount(g)
  e_all   <- igraph::ecount(g)
  g       <- igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  message(sprintf("  full network %d nodes / %d edges  ->  LCC %d nodes / %d edges",
                  n_all, e_all, igraph::vcount(g), igraph::ecount(g)))

  cm  <- igraph::cluster_louvain(g)
  mcc <- mcc_score(g)
  metrics <- c("degree", "betweenness", "closeness", "eigenvector",
               if (!is.null(mcc)) "MCC")

  tab <- tibble::tibble(
    gene             = igraph::V(g)$name,
    degree           = igraph::degree(g),
    degree_norm      = igraph::degree(g, normalized = TRUE),
    betweenness      = igraph::betweenness(g, normalized = TRUE),
    closeness        = igraph::closeness(g, normalized = TRUE),
    eigenvector      = igraph::eigen_centrality(g)$vector,
    clustering_coeff = igraph::transitivity(g, type = "local", isolates = "zero"),
    module           = as.integer(igraph::membership(cm))
  )
  if (!is.null(mcc)) tab$MCC <- mcc[tab$gene]

  tab <- tab |>
    dplyr::mutate(
      logFC     = lfc_map[gene],
      direction = dplyr::if_else(logFC > 0, "Up", "Down"),
      module_sz = as.integer(table(module)[as.character(module)]),
      dplyr::across(dplyr::all_of(metrics), minmax, .names = "{.col}_n")
    ) |>
    dplyr::mutate(
      hub_score      = rowMeans(dplyr::pick(dplyr::all_of(paste0(metrics, "_n")))),
      hub_components = paste(metrics, collapse = "+")
    ) |>
    dplyr::arrange(dplyr::desc(hub_score)) |>
    dplyr::mutate(
      hub_rank = dplyr::min_rank(dplyr::desc(hub_score)),
      is_hub   = hub_rank <= N_HUBS_REPORT
    )

  list(graph = g, edges = igraph::as_data_frame(g, "edges"), table = tab,
       communities = cm, mcc_used = !is.null(mcc),
       n_nodes_full = n_all, n_edges_full = e_all)
}

## ---- 4. Run per contrast ---------------------------------------------------

message("\ninteraction source: ", STRING_SOURCE,
        if (STRING_SOURCE == "local") "  (complete interactome, no query-size limit)" else "")

summaries <- list()

for (i in seq_len(nrow(CONTRASTS))) {
  cc <- CONTRASTS[i, ]
  sc <- paste0(cc$id, "_status")
  lc <- paste0(cc$id, "_logFC")

  degs    <- membership |> dplyr::filter(.data[[sc]] %in% c("Up", "Down"))
  lfc_map <- setNames(degs[[lc]], degs$Symbol)

  message(sprintf("\n%s (%s): %d DEGs", cc$id, cc$label, nrow(degs)))
  edges <- get_edges(degs$Symbol, cc$id)
  net   <- analyse_network(edges, lfc_map, cc$id)

  write_result(net$edges, file.path(DIR_PPI, sprintf("edges_%s.csv", cc$id)))
  write_result(net$table, file.path(DIR_PPI, sprintf("centrality_%s.csv", cc$id)))

  mod_summary <- net$table |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      n_nodes   = dplyr::n(),
      n_up      = sum(direction == "Up", na.rm = TRUE),
      pct_up    = round(100 * n_up / n_nodes, 1),
      top_nodes = paste(gene[order(-hub_score)][seq_len(min(6L, dplyr::n()))],
                        collapse = "; "),
      .groups   = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_nodes))
  write_result(mod_summary, file.path(DIR_PPI, sprintf("module_summary_%s.csv", cc$id)))

  n_up <- sum(net$table$direction == "Up", na.rm = TRUE)
  summaries[[cc$id]] <- tibble::tibble(
    contrast       = cc$id,
    label          = cc$label,
    n_degs         = nrow(degs),
    n_nodes        = igraph::vcount(net$graph),
    n_edges        = igraph::ecount(net$graph),
    n_nodes_before_lcc = net$n_nodes_full,
    n_edges_before_lcc = net$n_edges_full,
    n_modules      = length(unique(net$table$module)),
    n_up_nodes     = n_up,
    pct_up         = round(100 * n_up / igraph::vcount(net$graph), 1),
    modularity     = round(igraph::modularity(net$communities), 3),
    top_hub        = net$table$gene[1],
    hub_components = net$table$hub_components[1],
    mcc_exact      = net$mcc_used,
    string_source  = STRING_SOURCE,
    string_version = STRING_VERSION,
    score_min      = STRING_SCORE_MIN / 1000,
    retrieved_on   = as.character(Sys.Date())
  )

  message(sprintf("  %s: %d nodes, %d edges, %d modules, %.1f%% up, top hub = %s%s",
                  cc$id, igraph::vcount(net$graph), igraph::ecount(net$graph),
                  length(unique(net$table$module)),
                  100 * n_up / igraph::vcount(net$graph), net$table$gene[1],
                  if (!net$mcc_used) "  [MCC omitted]" else ""))
}

net_summary <- dplyr::bind_rows(summaries)
write_result(net_summary, file.path(DIR_PPI, "network_summary.csv"))

## Node and edge counts quoted anywhere downstream must come from this table:
## both are read from the same igraph object, after LCC extraction. The columns
## n_*_before_lcc are recorded so the effect of that step is auditable rather
## than inferred.
print(as.data.frame(
  dplyr::select(net_summary, contrast, n_degs, n_nodes, n_edges, n_modules,
                pct_up, top_hub, mcc_exact)
))

if (any(!net_summary$mcc_exact)) {
  message("\nNOTE: MCC was omitted for ",
          paste(net_summary$contrast[!net_summary$mcc_exact], collapse = ", "),
          ". Their composite hub scores average four metrics, not five; the ",
          "hub_components column records this per row. Do not describe the score ",
          "as five-metric for those contrasts.")
}

message("\n05_ppi_networks.R complete.")
