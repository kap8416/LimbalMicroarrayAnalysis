## =============================================================================
## run_all.R
## Sequential driver for the full pipeline. Run from the repository root:
##   Rscript R/run_all.R
##
## Steps 01 and 05 require network access (GEO and the STRING REST API).
## Set OFFLINE <- TRUE in R/00_config.R to reuse cached downloads in data/raw/.
## =============================================================================

steps <- c(
  "R/01_download_data.R",
  "R/02_integration.R",
  "R/03_differential_expression.R",
  "R/04_functional_enrichment.R",
  "R/05_ppi_networks.R",
  "R/00b_build_tf_census.R",
  "R/06_mra_rewiring.R",
  "R/07_core_signature_null.R"
)

t0 <- Sys.time()

for (s in steps) {
  message("\n", strrep("=", 78))
  message("RUNNING  ", s)
  message(strrep("=", 78))
  t1 <- Sys.time()
  source(s, echo = FALSE)
  message(sprintf("--- %s finished in %.1f min", s,
                  as.numeric(difftime(Sys.time(), t1, units = "mins"))))
}

message("\n", strrep("=", 78))
message(sprintf("PIPELINE COMPLETE in %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
message(strrep("=", 78))

writeLines(capture.output(sessionInfo()), "env/sessionInfo.txt")
message("session information written to env/sessionInfo.txt")
message("\nNext: generate figures with")
message("  for f in figures/scripts/fig*.py; do python \"$f\"; done")
