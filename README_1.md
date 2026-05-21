# LimbalMicroarrayAnalysis

Reproducible analysis pipeline for the integrative transcriptomic characterization of human limbal epithelial identity.

> **Aviña-Padilla et al.** — *Integrative transcriptomic analysis of the human limbal epithelium reveals a conserved niche-dependent regulatory architecture* — manuscript in preparation.

---

## Overview

This repository contains all R scripts, processed data, and output figures associated with the systems-level analysis of limbal epithelial identity using two public human microarray datasets (GSE38190 / GPL570 and GSE56421 / GPL6244). The pipeline compares limbal epithelium against corneal epithelium, conjunctival epithelium, and cultured limbal epithelial cells across three limbus-centered contrasts.

```
GEO raw data (GSE38190 + GSE56421)
        │
        ▼
Batch correction (quantile norm → cornea-anchored shift → ComBat)
        │
        ▼
Gene-level collapse (MAD-based isoform selection, 24,002 genes)
        │
        ├──▶ Differential Expression (limma)
        │         └──▶ Core limbal signature (32 genes)
        │
        ├──▶ Functional Enrichment (ORA — GO BP / KEGG)
        │
        ├──▶ PPI Network (STRING v12.0 / Cytoscape)
        │
        └──▶ Co-expression & MRA (corto)
                  └──▶ Regulatory Rewiring
```

---

## Repository Structure

```
LimbalMicroarrayAnalysis/
│
├── scripts/
│   ├── DEA_limbal_v3_FINAL_repro.R          # Differential expression analysis
│   ├── functional_enrichment_repro.R         # ORA — GO BP and KEGG
│   ├── corto_mra_rewiring_repro.R            # Co-expression network, MRA, rewiring
│   └── ppi_network.R                         # PPI network construction (STRING)
│
├── data/
│   ├── combined_expression_matrix_corrected_v2.csv   # Batch-corrected matrix (45,378 transcripts)
│   ├── combined_expression_matrix_gene_level.csv     # Gene-level matrix (24,002 genes)
│   ├── probe_to_gene_mapping.csv                     # RefSeq → HGNC symbol map (MAD-based)
│   └── sample_metadata.csv                           # Sample annotations
│
├── results/
│   ├── DEA/
│   │   ├── DEA_Limbus_vs_Cornea.csv
│   │   ├── DEA_Limbus_vs_Cornea_UP.csv
│   │   ├── DEA_Limbus_vs_Cornea_DOWN.csv
│   │   ├── DEA_Limbus_vs_Conjunctiva.csv
│   │   ├── DEA_Limbus_vs_Conjunctiva_UP.csv
│   │   ├── DEA_Limbus_vs_Conjunctiva_DOWN.csv
│   │   ├── DEA_Limbus_vs_Cultured.csv
│   │   ├── DEA_Limbus_vs_Cultured_UP.csv
│   │   ├── DEA_Limbus_vs_Cultured_DOWN.csv
│   │   ├── DEA_counts_summary.csv
│   │   └── Limbal_core_signature_overlaps.csv
│   │
│   ├── enrichment/
│   │   └── enrichment_summary_all.csv
│   │
│   ├── corto/
│   │   ├── regulon.rds
│   │   ├── regulon_table.csv
│   │   ├── MRA_Limbus_vs_Cornea.csv
│   │   ├── MRA_Limbus_vs_Conjunctiva.csv
│   │   ├── MRA_Limbus_vs_Cultured.csv
│   │   └── rewiring_summary.csv
│   │
│   └── figures/
│       ├── PCA_final_epithelial_v3.png
│       ├── Volcano_Limbus_vs_Cornea.png
│       ├── Volcano_Limbus_vs_Conjunctiva.png
│       ├── Volcano_Limbus_vs_Cultured.png
│       ├── Heatmap_Limbus_vs_Cornea.pdf
│       ├── Heatmap_Limbus_vs_Conjunctiva.pdf
│       ├── Heatmap_Limbus_vs_Cultured.pdf
│       ├── UpSet_limbal_upregulated.pdf
│       ├── regulon_network.png
│       └── rewiring_network.png
│
└── session/
    └── sessionInfo_DEA_v3.txt
```

---

## Data Sources

| GEO Accession | Platform | Samples used |
|---------------|----------|-------------|
| GSE38190 | GPL570 (Affymetrix HG-U133 Plus 2.0) | Limbus n=4, Cornea n=3, Conjunctiva n=3 (GSM724094 excluded) |
| GSE56421 | GPL6244 (Affymetrix HuGene-1_0-st) | Cornea n=3, Cultured limbal epithelial cells n=3 (fibroblasts excluded) |

Raw data are publicly available from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

---

## Reproducibility

All scripts include deterministic fixes to ensure identical results across runs:

- `set.seed(42)` + `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` at pipeline start
- MAD-based isoform collapse with lexicographic tiebreaker (`order(-mad, refseq)`)
- Deterministic annotation map construction (`sort(names())` before deduplication)
- MD5 checksum written to `reproducibility_checksum.txt` on each run

To verify reproducibility, run `DEA_limbal_v3_FINAL_repro.R` twice and confirm identical values in `reproducibility_checksum.txt`.

---

## Execution Order

Scripts must be run in the following order:

```r
# 1. Differential expression + gene-level matrix export
source("scripts/DEA_limbal_v3_FINAL_repro.R")

# 2. Functional enrichment
source("scripts/functional_enrichment_repro.R")

# 3. Co-expression network, MRA, and regulatory rewiring
source("scripts/corto_mra_rewiring_repro.R")

# 4. PPI network (requires STRING output loaded in Cytoscape)
source("scripts/ppi_network.R")
```

> **Note:** `corto_mra_rewiring_repro.R` reads `combined_expression_matrix_gene_level.csv` exported by step 1. Always run step 1 first or ensure this file is present in the working directory.

---

## Dependencies

All analyses were performed in R. Required packages:

| Package | Version | Source | Purpose |
|---------|---------|--------|---------|
| limma | ≥3.54 | Bioconductor | Differential expression |
| sva | ≥3.46 | Bioconductor | ComBat batch correction |
| corto | ≥1.1 | CRAN | Co-expression network + MRA |
| clusterProfiler | ≥4.6 | Bioconductor | Functional enrichment |
| org.Hs.eg.db | ≥3.16 | Bioconductor | Gene annotation |
| hgu133plus2.db | ≥3.2 | Bioconductor | GPL570 annotation |
| hugene10sttranscriptcluster.db | ≥8.8 | Bioconductor | GPL6244 annotation |
| ggplot2 | ≥3.4 | CRAN | Visualization |
| ggrepel | ≥0.9 | CRAN | Volcano plot labels |
| pheatmap | ≥1.0 | CRAN | Heatmaps |
| UpSetR | ≥1.4 | CRAN | UpSet plots |
| dplyr | ≥1.1 | CRAN | Data manipulation |
| digest | ≥0.6 | CRAN | Reproducibility checksums |

Install all Bioconductor dependencies:

```r
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "limma", "sva", "clusterProfiler",
  "org.Hs.eg.db", "hgu133plus2.db",
  "hugene10sttranscriptcluster.db", "AnnotationDbi"
))
```

See `session/sessionInfo_DEA_v3.txt` for the exact R and package versions used to generate the published results.

---

## Key Results

| Contrast | Upregulated | Downregulated |
|----------|-------------|---------------|
| Limbus vs Cornea | 1,865 | 934 |
| Limbus vs Conjunctiva | 112 | 36 |
| Limbus vs Cultured | 1,009 | 548 |

**Core limbal signature (32 genes, conserved across all three contrasts):**
ABI3BP, ADARB1, ANGPT1, BMP4, C14orf132, CCDC102B, CDH11, CEP126, COL8A1, COL8A2, CYTL1, EYA4, FGF18, FIBIN, FZD7, GRP, HACD1, ITGA10, ITGBL1, JCAD, KCNE4, LRRC3B, NTRK3, PLPP3, POU3F3, PRR16, PTGFR, SHC3, SHOX, SOD3, THBS4, ZEB1

**Dominant limbal master regulators:** ZEB2, ETS1, ZEB1, MEF2C

**Top rewired TFs:** MEF2C, MITF, ZEB2, IKZF1, FLI1, NFIX, SMAD5

---

## Citation

> Aviña-Padilla et al. *Integrative transcriptomic analysis of the human limbal epithelium reveals a conserved niche-dependent regulatory architecture.* Manuscript in preparation.

GEO accessions: [GSE38190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38190) · [GSE56421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56421)

---

## License

Scripts are released under the [MIT License](LICENSE). Processed data are shared for reproducibility purposes under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
