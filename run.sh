#!/usr/bin/env bash
# =============================================================================
# run.sh — one-command runner for the whole pipeline
#
# USAGE (from anywhere):
#     bash "/full/path/to/LimbalMicroarrayAnalysis/run.sh"
#
# Or, from inside the repository:
#     ./run.sh
#
# Do NOT paste the contents of the R scripts into an interactive R console. They
# depend on being sourced in order from the repository root, so a paste loses both
# the working directory and the objects defined by 00_config.R. This script sets
# the working directory itself, so it cannot go wrong.
#
# Options:
#     ./run.sh --from-existing-matrix   use data/raw/combined_expression_matrix_corrected_v2.csv
#                                       instead of regenerating it from GEO (recommended
#                                       for the first clean re-run; see docs/AUDIT.md 3.7)
#     ./run.sh --skip-string-download   assume the STRING flat files are already present
#     ./run.sh --figures-only           regenerate figures from existing results/
# =============================================================================

set -uo pipefail

# Always operate from the directory containing this script.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO" || exit 1

USE_EXISTING=0
SKIP_STRING=0
FIGURES_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --from-existing-matrix) USE_EXISTING=1 ;;
    --skip-string-download) SKIP_STRING=1 ;;
    --figures-only)         FIGURES_ONLY=1 ;;
    -h|--help)              sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "Unknown option: $arg  (try --help)"; exit 2 ;;
  esac
done

mkdir -p logs env
LOG="logs/run_$(date +%Y%m%d_%H%M%S).log"

# Everything below is written to the terminal AND to $LOG.
exec > >(tee -a "$LOG") 2>&1

rule() { printf '%s\n' "------------------------------------------------------------------------"; }

echo "LimbalMicroarrayAnalysis — pipeline run"
rule
echo "repository : $REPO"
echo "log        : $LOG"
echo "started    : $(date)"
echo "host       : $(uname -srm)"
rule

# ---------------------------------------------------------------------------
# 1. Prerequisites
# ---------------------------------------------------------------------------
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not on PATH."; exit 1; }

# Verify the pipeline itself is present before doing any downloading or installing.
# The scripts live in R/ next to this file; if the repository was copied or moved
# partially, fail here with a useful message rather than midway through a run.
MISSING=()
for f in 00_config.R 01_download_data.R 03_differential_expression.R \
         04_functional_enrichment.R 05_ppi_networks.R 06_mra_rewiring.R \
         07_core_signature_null.R; do
  [ -s "$REPO/R/$f" ] || MISSING+=("R/$f")
done
if [ "${#MISSING[@]}" -gt 0 ]; then
  echo "ERROR: the pipeline scripts are missing from $REPO/R/"
  printf '  %s\n' "${MISSING[@]}"
  echo
  echo "Every script must sit in the R/ directory beside run.sh. If you saved them"
  echo "individually (for example into ~/Downloads), move them back:"
  echo "    mkdir -p \"$REPO/R\""
  echo "    mv ~/Downloads/{00,01,02,02b,03,04,05,06,07,99}*.R \"$REPO/R/\""
  echo "    mv ~/Downloads/run_all.R \"$REPO/R/\""
  exit 1
fi
echo "pipeline scripts: all present in R/"
Rscript -e 'cat(R.version.string, "\n")'

PY=""
for cand in python3 python; do
  command -v "$cand" >/dev/null 2>&1 && { PY="$cand"; break; }
done
[ -n "$PY" ] && "$PY" --version || echo "WARNING: no python found; figures will be skipped."

echo
echo "Checking R packages (installing anything missing) ..."
Rscript - <<'RS'
cran <- c("tidyverse", "igraph", "httr", "jsonlite", "purrr", "tibble", "readr", "dplyr")
bioc <- c("GEOquery", "limma", "sva", "clusterProfiler", "org.Hs.eg.db",
          "AnnotationDbi", "Biobase", "corto")

need <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) {
  cat("installing from CRAN:", paste(need, collapse = ", "), "\n")
  install.packages(need, repos = "https://cloud.r-project.org", quiet = TRUE)
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)

need <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) {
  cat("installing from Bioconductor:", paste(need, collapse = ", "), "\n")
  BiocManager::install(need, ask = FALSE, update = FALSE)
}

missing <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("\nERROR: still missing:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
}
cat("all R packages available\n")
RS
[ $? -ne 0 ] && { echo "ERROR: R package setup failed. See $LOG"; exit 1; }

# ---------------------------------------------------------------------------
# 2. STRING flat files
# ---------------------------------------------------------------------------
if [ "$FIGURES_ONLY" -eq 0 ] && [ "$SKIP_STRING" -eq 0 ]; then
  echo
  echo "Checking STRING v12.0 flat files ..."
  mkdir -p data/raw
  for f in 9606.protein.links.v12.0.txt.gz 9606.protein.aliases.v12.0.txt.gz; do
    if [ -s "data/raw/$f" ]; then
      echo "  present: $f ($(du -h "data/raw/$f" | cut -f1))"
    else
      case "$f" in
        *links*)   sub="protein.links.v12.0" ;;
        *aliases*) sub="protein.aliases.v12.0" ;;
      esac
      echo "  downloading $f ..."
      curl -fL --retry 3 -o "data/raw/$f" \
        "https://stringdb-downloads.org/download/$sub/$f" \
        || { echo "  ERROR: download failed. Fetch it manually into data/raw/ and re-run with --skip-string-download"; exit 1; }
    fi
  done
fi

# ---------------------------------------------------------------------------
# 3. Clear stale outputs
# ---------------------------------------------------------------------------
# Result files ingested from an earlier analysis, or left by a partial run, can be
# read-only or locked and would then block the new run midway. They are also a
# provenance hazard: a directory holding a mix of old and new tables cannot be
# reasoned about. The originals remain in data/raw/legacy/ and provenance/.
if [ "$FIGURES_ONLY" -eq 0 ]; then
  echo
  echo "Clearing stale result files ..."
  n=0
  for d in results/de results/enrichment results/ppi results/mra data/processed; do
    [ -d "$d" ] || continue
    for f in "$d"/*; do
      [ -f "$f" ] || continue
      case "$(basename "$f")" in .gitkeep) continue ;; esac
      chmod u+w "$f" 2>/dev/null
      command -v chflags >/dev/null 2>&1 && chflags nouchg "$f" 2>/dev/null
      rm -f "$f" 2>/dev/null && n=$((n+1)) || echo "  WARNING: could not remove $f"
    done
  done
  echo "  removed $n file(s)"
fi

# ---------------------------------------------------------------------------
# 4. Pipeline
# ---------------------------------------------------------------------------
if [ "$FIGURES_ONLY" -eq 0 ]; then
  if [ "$USE_EXISTING" -eq 1 ]; then
    MATRIX="data/raw/combined_expression_matrix_corrected_v2.csv"
    [ -s "$MATRIX" ] || { echo "ERROR: $MATRIX not found. Copy it there first."; exit 1; }
    STEPS=(R/01_download_data.R R/02b_use_existing_matrix.R)
    echo
    echo "Integration: using the existing harmonized matrix (02b)."
    echo "  Step 01 still queries GEO to verify the sample accessions."
  else
    STEPS=(R/01_download_data.R R/02_integration.R)
    echo
    echo "Integration: regenerating the matrix from GEO (02)."
  fi
  STEPS+=(R/03_differential_expression.R
          R/04_functional_enrichment.R
          R/05_ppi_networks.R
          R/00b_build_tf_census.R
          R/06_mra_rewiring.R
          R/07_core_signature_null.R)

  for s in "${STEPS[@]}"; do
    echo
    rule; echo "RUNNING  $s"; rule
    t0=$(date +%s)
    if Rscript "$REPO/$s"; then
      echo "--- OK  $s  ($(( ($(date +%s) - t0) / 60 ))m $(( ($(date +%s) - t0) % 60 ))s)"
    else
      echo
      echo "########################################################################"
      echo "FAILED at $s"
      echo "Everything above is in: $LOG"
      echo "Send that log and we can fix it. Earlier steps completed and their"
      echo "outputs are in results/, so the run can resume from here."
      echo "########################################################################"
      exit 1
    fi
  done
fi

# ---------------------------------------------------------------------------
# 5. Figures
# ---------------------------------------------------------------------------
if [ -n "$PY" ]; then
  echo
  rule; echo "FIGURES"; rule
  "$PY" -m pip install -q -r env/requirements.txt 2>/dev/null
  ( cd figures/scripts || exit 1
    "$PY" reconstruct_c2_edges.py 2>/dev/null || true
    for f in fig*.py; do
      if "$PY" "$f" >/dev/null 2>&1; then echo "  OK       $f"
      else                                echo "  skipped  $f  (missing input — run it directly for the reason)"; fi
    done )
fi

# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------
echo
rule; echo "SUMMARY"; rule
for f in results/de/integration_qc.csv results/de/DEA_counts_summary.csv \
         results/de/duplicate_profile_report.csv results/ppi/network_summary.csv \
         results/de/core_signature_null.txt; do
  [ -s "$f" ] && { echo; echo "== $f"; cat "$f"; }
done

echo
echo "figures written:"
ls -1 figures/output/*.pdf 2>/dev/null | sed 's|^|  |' || echo "  none"

Rscript -e 'writeLines(capture.output(sessionInfo()), "env/sessionInfo.txt")' 2>/dev/null \
  && echo && echo "session info: env/sessionInfo.txt"

echo
echo "finished: $(date)"
echo "full log: $LOG"
