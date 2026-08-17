<h1 align="center">LimbalMicroarrayAnalysis</h1>

<p align="center">
  <strong>Comparative systems-level transcriptomics of human limbal epithelial identity</strong><br>
  <em>A reproducible, auditable reanalysis of two public ocular-surface microarray datasets</em>

---

## What this is

The corneal epithelium depends on the limbus for lifelong renewal, yet the
regulatory architecture that defines limbal identity is usually described one
marker at a time. This repository reanalyses two public human ocular-surface
microarray datasets as a **single harmonized matrix**, and interrogates limbal
identity through five analytical layers rather than a differential-expression list
alone: differential expression, functional over-representation, protein–protein
interaction topology, master regulator inference, and condition-specific
co-expression rewiring.

The central result is that limbal identity behaves as a **niche-dependent
regulatory attractor**: a self-reinforcing network state anchored on FN1-centred
extracellular-matrix signalling, which is destabilized once cells leave their
native niche.


---

## Design

Three contrasts, all anchored on the limbus, so a positive log₂ fold change always
means higher expression in limbus:

| | Comparison | Platforms | \|log₂FC\| | Interpretation |
|---|---|---|---|---|
| **C1** | limbus vs corneal epithelium | cross-platform | ≥ 1.0 | progenitor state → terminal differentiation |
| **C2** | limbus vs conjunctival epithelium | same platform (GPL570) | ≥ 0.5 | adjacent but distinct lineages — **most reliable contrast** |
| **C3** | limbus vs cultured limbal epithelial cells | confounded — **exploratory** | ≥ 1.0 | consequences of leaving the niche |

All contrasts use FDR ≤ 0.05 (Benjamini–Hochberg). C3 is fully confounded with
platform of origin (limbus = GPL570, cultured cells = GPL6244) and is reported as
exploratory throughout; it requires independent experimental validation.

### Data

Raw data are **not redistributed**; `R/01_download_data.R` retrieves them from GEO.

| Accession | Platform | Retained groups |
|---|---|---|
| [GSE38190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38190) | GPL570 · HG-U133 Plus 2.0 | limbus 4 · cornea 3 · conjunctiva 3 |
| [GSE56421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56421) | GPL6244 · HuGene-1_0-st | cornea 3 · cultured limbal epithelial cells 3 |

**16 samples across 4 groups** after curation. Twenty samples are available; four
are excluded and the exclusions are enforced by assertion in every layer:

- `GSM1361191–93` — limbal fibroblasts, non-epithelial.
- `GSM724094` — conjunctival outlier (within-group *r* = 0.92 versus 0.98–0.99).

Cornea is the only group present on both platforms and therefore the sole anchor
for cross-platform harmonization.

---

## Quickstart

```bash
git clone https://github.com/kap8416/LimbalMicroarrayAnalysis.git
cd LimbalMicroarrayAnalysis

# R dependencies
Rscript -e 'install.packages("renv"); renv::restore()'

# Python dependencies (figures only)
python -m pip install -r env/requirements.txt

# One-time: the complete STRING interactome (~90 MB, see note below)
curl -L -o data/raw/9606.protein.links.v12.0.txt.gz \
  https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
curl -L -o data/raw/9606.protein.aliases.v12.0.txt.gz \
  https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz

# Full pipeline
Rscript R/run_all.R

# Figures
for f in figures/scripts/fig*.py; do python "$f"; done
```

Step 01 queries GEO. Set `OFFLINE <- TRUE` in `R/00_config.R` to reuse cached
downloads in `data/raw/`.

> [!NOTE]
> **Why the local STRING files rather than the API.** The REST `/network` endpoint
> returns only interactions *among the identifiers submitted in that request*. C1
> has 2,799 differentially expressed genes, above the practical single-request
> limit, so any API-based approach must either chunk the list — which silently
> discards every interaction spanning two chunks — or fail. The local flat files
> contain the complete interactome, have no size limit, and are versioned, so the
> query is exactly reproducible. `R/05_ppi_networks.R` refuses to chunk and stops
> with instructions instead.

---

## Pipeline

```
R/
├── 00_config.R                   paths · thresholds · contrasts · palettes · QC gates
├── 01_download_data.R            GEO retrieval · sample curation · GSM assertions
├── 02_integration.R              quantile norm → cornea-anchored shift → ComBat → MAD collapse
├── 03_differential_expression.R  limma · 3 contrasts · core signature
├── 04_functional_enrichment.R    clusterProfiler GO-BP + KEGG, matrix-restricted universe
├── 05_ppi_networks.R             STRING v12.0 · igraph · Louvain · exact MCC
├── 06_mra_rewiring.R             corto regulons · MRA · |Δρ| rewiring
├── 07_core_signature_null.R      permutation null for the core signature
├── 99_ingest_legacy_results.R    map pre-existing result files onto the canonical names
└── run_all.R                     sequential driver, writes env/sessionInfo.txt
```

Each script sources `00_config.R`, so thresholds and paths are defined once. Every
step writes tidy CSVs to `results/`, which the Python figure scripts consume; no
figure reads an intermediate R object.

### Method summary

**Integration.** Probe signals are summarized as RefSeq-level estimates and
collapsed to gene level by retaining, per symbol, the isoform with the highest
median absolute deviation. Because cornea is the only shared group, harmonization
proceeds in three steps: within-platform quantile normalization
(`limma::normalizeBetweenArrays`), cornea-anchored cross-platform shift correction,
and residual `ComBat` (`sva`, `mod = ~is_Cornea` — a biological covariate *is*
protected, which matters because group and batch are partially confounded).
Gated on cross-platform corneal concordance ≥ 0.95.

**Differential expression.** `limma`, no-intercept design, `eBayes(trend = TRUE)`.
Not voom: voom models the mean–variance relationship of counts and does not apply
to normalized log₂ array intensities.

**Enrichment.** `clusterProfiler` `enrichGO` (GO-BP) and `enrichKEGG`, up- and
downregulated sets separately. The universe is restricted to genes present in the
batch-corrected matrix — using the whole genome would inflate significance, since
the measurable space is bounded by the platform intersection.

**PPI.** STRING v12.0, combined confidence ≥ 0.70, functional network. Topology in
`igraph`; hubs ranked by the mean of five min–max normalized centralities
(degree, betweenness, closeness, eigenvector, **exact** maximum clique centrality);
largest connected component retained; communities by Louvain.

**Master regulators.** `corto` on the top 5,000 most variable genes with annotated
TFs as centroids, 1,000 bootstraps, *p* ≤ 1 × 10⁻⁶, `minsize = 10`. NES scores
**regulon activity, not transcript abundance** — a distinction the figures state
explicitly, because PAX6 is constitutively expressed while its regulon activity
falls.

**Rewiring.** Per-group Spearman co-expression; edges with |Δρ| ≥ 0.6 between
conditions classified as rewired.

---

## Figures

Publication-ready PDF (vector, TrueType-embedded) plus 600 dpi PNG in
`figures/output/`. Style and palette are centralized in `figures/scripts/fig_style.py`
and mirrored in `R/00_config.R`.

| | Content | Status |
|---|---|---|
| Fig. 1 | Analytical pipeline | ✅ |
| Fig. 2 | DEG burden + volcano panels | ✅ |
| Fig. 3 | Set overlap + 32-gene core signature | ✅ |
| Fig. 4 | GO-BP / KEGG dot plots | ⏳ awaiting enrichment tables |
| Fig. 5 | C1 interactome, FN1 hub, module landscape | ⚠️ regenerate after the PPI fix |
| Fig. 6 | C2 circuit + C3 interactome | ⚠️ regenerate after the PPI fix |
| Fig. 7 | Regulatory attractor model | ✅ |
| Fig. 8 | Master regulators + rewiring | ⏳ awaiting MRA tables |
| Suppl. 1 | PCA + sample correlation QC | ⏳ awaiting the harmonized matrix |

Scripts awaiting input exit with a message naming the missing file and the R step
that produces it; they have been tested against fixtures of the expected schema.
Legends are in [`docs/FIGURES.md`](docs/FIGURES.md).

---

## Reproducibility status

An independent audit of the original analysis is recorded in
[`docs/AUDIT.md`](docs/AUDIT.md). It separates checks that **passed** from design
concerns and from code defects, and states explicitly what could not be verified.

**Verified sound.** DEG calls match the declared thresholds exactly (0 discrepancies
in 3,448 genes × 3 contrasts); the composite hub score equals the mean of its
normalized components to 10⁻¹⁶; Benjamini–Hochberg adjustment reproduces
independently; the C3 *p*-value histogram is well behaved, with no signature of
over-correction. The published DEG counts, the 32-gene core signature, and the
FN1/IL6/LAMB3 hub values all reproduce from the deposited tables.

**Corrected here.**

| # | Defect | Effect | Fix |
|---|---|---|---|
| [3.3](docs/AUDIT.md) | STRING query chunked at 1,500 genes; cross-chunk interactions never queried | C1 and C3 networks missing **12%** and **5%** of edges; all their topology metrics computed on incomplete graphs | local interactome; chunking now refused |
| [3.4](docs/AUDIT.md) | `MCC` replaced by `degree²` when *n* > 500 | two of five composite terms were functions of degree, weighting it 2/5 not 1/5 | exact MCC, or the term is dropped and the omission recorded |
| [3.5](docs/AUDIT.md) | MRA and rewiring run on 17 samples, retaining the excluded outlier | conjunctiva *n* = 4 in the regulatory layer, *n* = 3 in the DE layer | sample-set assertions in every layer |
| [2.1](docs/ERRATA.md) | C2 top-DEG values are C1 fold changes (7 genes; CRTAC1 and GJB6 sign-reversed) | manuscript text | corrected values listed |
| [2.2](docs/ERRATA.md) | C2 network reported as 10 edges | topology is 4 edges, derived and validated from the centrality vectors | `figures/scripts/reconstruct_c2_edges.py` |

**Open.** The script that produced the harmonized matrix
(`combined_expression_matrix_corrected_v2.csv`) has not been recovered, so the
corneal concordance of *r* = 0.953, the PCA structure and the isoform collapse
remain unverified against the original run. `R/02_integration.R` implements the step
as the Methods describe it and can regenerate the matrix.

Also open, and worth addressing in revision rather than defending: the three
contrasts **share the same four limbal samples**, so they are not independent tests
and the core signature's "stability across tissue contexts" is stronger than the
design supports. `R/07_core_signature_null.R` provides a within-platform
permutation null that respects that dependency.

`provenance/` archives the original scripts unmodified, with a timestamp chain
linking each to the tables it produced.

---

## Repository layout

```
├── R/                  analysis pipeline
├── figures/
│   ├── scripts/        figure generation (matplotlib)
│   └── output/         PDF + PNG, 600 dpi
├── results/
│   ├── de/             limma tables · core signature · counts
│   ├── enrichment/     GO-BP and KEGG
│   ├── ppi/            edges · centralities · modules · network summary
│   └── mra/            regulons · NES · rewiring
├── data/
│   ├── raw/            GEO downloads + STRING flat files (gitignored)
│   └── processed/      harmonized matrices (gitignored)
├── provenance/         original scripts, archived for the record
├── docs/
│   ├── FIGURES.md      legends · panel inventory · generation status
│   ├── ERRATA.md       manuscript values checked against the tables
│   └── AUDIT.md        soundness of the analyses
└── env/                renv.lock · requirements.txt · sessionInfo.txt
```

---

## Key quantities

Re-derived from the deposited tables. Values that did not reconcile with the
manuscript draft are marked and explained in `docs/ERRATA.md`.

| | |
|---|---|
| Gene-level features after collapse | 24,002 |
| Cross-platform corneal concordance | *r* = 0.953 |
| DEGs — C1 / C2 / C3 | 1,865↑ 934↓ · 112↑ 36↓ · 1,009↑ 548↓ |
| Core limbal signature | 32 named genes (of 49 significant features) |
| Dominant C1 hub | FN1 · degree 90 · log₂FC 4.81 |
| Most central C2 node | LAMB3 · betweenness 0.833 |
| C2 network | 5 nodes, **4** edges (not 10 — `ERRATA.md` §2.2) |
| Top C1 master regulator | ETS1 · NES 7.27 |
| Most rewired regulator | MEF2C · 485 edges |
| Genes down in all three contrasts | **1** (versus 49 up — see `AUDIT.md` §2.6) |

---

## Citation

If you use this code, please cite both the article and the repository. Machine-readable
metadata is in [`CITATION.cff`](CITATION.cff).

> Aviña-Padilla K, Cárdenas-López LE, Valencia-Lozano E, Castro-Muñozledo F (2026).
> *Human limbal epithelial identity is maintained by a niche-dependent regulatory
> attractor anchored on FN1.* (submitted).

---

## Data and code availability

All primary data are public: [GSE38190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38190)
and [GSE56421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56421).
Interaction data are STRING v12.0. Derived tables are in `results/`. Code in this
repository is MIT-licensed; the GEO and STRING data remain under their original
terms.

## Contact

| | |
|---|---|
| **Katia Aviña-Padilla** | Departamento de Biología Celular, CINVESTAV-IPN |
| **Federico Castro-Muñozledo** (corresponding) | federico.castro@cinvestav.mx |

## License

MIT — see [`LICENSE`](LICENSE).
