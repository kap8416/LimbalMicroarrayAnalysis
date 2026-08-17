# Environment

## R (analysis pipeline)

R >= 4.3 with Bioconductor >= 3.18.

```r
install.packages("renv")
renv::restore()          # from env/renv.lock
```

`renv.lock` is generated on the analysis machine, not committed by hand:

```r
renv::init()             # once, at the repository root
renv::snapshot(lockfile = "env/renv.lock")
```

Packages required by the pipeline:

| Step | Packages |
|---|---|
| 01 | GEOquery, tidyverse |
| 02 | limma, sva, Biobase |
| 03 | limma, purrr |
| 04 | clusterProfiler, org.Hs.eg.db, AnnotationDbi |
| 05 | httr, jsonlite, igraph |
| 06 | corto, org.Hs.eg.db, parallel |

## Python (figures)

Python >= 3.10.

```bash
python -m pip install -r env/requirements.txt
```

No plotting libraries beyond matplotlib are required: the network layout is
implemented directly in `figures/scripts/ppi_common.py` so that networkx is not a
dependency.

## Recording the session

`R/run_all.R` writes `env/sessionInfo.txt` at the end of a successful run. Commit
that file alongside any change to the results.
