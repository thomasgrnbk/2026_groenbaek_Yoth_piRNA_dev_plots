# Single TE Insertion Cut&Tag Heatmaps

Scripts for visualizing Cut&Tag (CUT&Tag) histone modification signal across individual transposable element (TE) insertions in *Drosophila melanogaster*, related to Figure 2D (and Supplementary Figure 9) in the associated paper.

**Publication:** *[link to be added upon publication]*  
**Raw sequencing data:** GEO accession *[to be added upon publication]*  
**Processed data:** [FTP server](http://ftp.genome-ftp.mbg.au.dk/public/PAN/Data_piDev_paper/cutntag_single_insertions_processed_data/) (public access)

---

## Biological context

This analysis examines how piRNA pathway activity (via *piwi* knockdown) affects histone modification states at individual TE insertions in the *Drosophila* germline. Chromatin state is assessed by Cut&Tag profiling of three histone marks:

| Mark | Biological association |
|------|----------------------|
| H3K9me3 | Heterochromatin / TE silencing |
| H3K4me3 | Active promoters / transcription |
| H3K27me3 | Polycomb repression (supplemental) |

Two genetic backgrounds / cell types are compared:

| Label | Background | Cell type | Knockdown driver |
|-------|-----------|-----------|-----------------|
| **BamMut** | *bam* mutant | GSC-like cells (sorted from ovary cysts) | *nanos*-Gal4 (early germline) |
| **MTD** / fullOvary | Full ovary | 512c-stage nurse cell nuclei (sorted) | MTD-Gal4 (germline development) |

Each condition is compared between:
- **control** – *white* shRNA (white KD)
- **knockdown** – *piwi* shRNA (piwi KD)

Signal values represent the **average CPM+1 across 2 biological replicates**.

---

## TE annotation

TEs were identified from two complementary sources:

- **Reference TEs** – identified in an in-house nanopore-sequenced genome assembly using [EDTA v2.2](https://github.com/oushujun/EDTA) with a *Drosophila* consensus TE library as reference.
- **Non-reference TEs** – additional insertions found by running [ngs_te_mapper2 v1.0.2](https://github.com/bergmanlab/ngs_te_mapper2) on all Cut&Tag FASTQ files from the fly crosses.

The combined catalog contains **1,467 single TE insertions** spanning the *D. melanogaster* genome.

---

## TE families and groupings

When plotting subsets (`SUB_TEs = TRUE`), TEs are grouped by their germline expression profile:

| Group | TE families | Color |
|-------|------------|-------|
| Somatic | 17.6, 412, Tirant | Gray |
| Germline early | mdg3, copia, accord | Orange |
| Germline late | 3S18, Burdock, blood | Blue |

---

## Filtering pipeline

The following filters are applied before heatmap generation (configurable at the top of `single_TE_insertion_hm.R`):

### 1. Euchromatin filter (`EUCH_FILT`)
Retains only TE insertions in euchromatic regions by requiring a minimum distance to the nearest neighboring TE insertion.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_NEXT_TE` | 1500 nt | Minimum distance to the nearest neighboring TE |
| `MIN_TE_COV` | 0 | Minimum TE consensus coverage (reference TEs only) |

### 2. H3K4me3 filter (`APPLY_H3K4_FILTER`)
Removes TE insertions that carry constitutive H3K4me3 signal unrelated to piRNA pathway activity — i.e., insertions where H3K4me3 is already high in control and does not change substantially upon *piwi* knockdown.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `H3K4_THRESHOLD` | 2 | Max allowed average upstream H3K4me3 CPM+1 in control |
| `H3K4_MAX_FC` | 2 | Minimum fold-change (piwi KD / control) to retain insertion |
| `FILTER_MODE` | `"union"` | `"union"`: filter if fails in ANY dataset; `"intersection"`: filter only if fails in ALL datasets |

Filtering is applied on raw CPM+1 values **before** log2 transformation to ensure unbiased fold-change estimates.

### Filtering summary (last run)
- Starting insertions (post-euchromatin filter): **420**
- Filtered by H3K4me3 (union of MTD + BamMut): **30**
- Final insertions in heatmap: **390** (281 reference + 109 non-reference)

---

## Heatmap details

- **Window:** ±500 bp around the TE insertion site
- **Bin size:** 50 bp (20 bins upstream, 20 bins downstream)
- **Signal transformation:** log2(CPM+1) for heatmap display
- **Row ordering:** by TE group → TE family → H3K9me3 signal (descending) in control
- **LFC heatmaps:** log2(piwi KD / white KD), computed from pre-log2 matrices

Color scales:

| Track | Color range | Scale |
|-------|------------|-------|
| H3K9me3 | white → red | 0–6 log2(CPM+1) |
| H3K4me3 | white → blue | 0–6 log2(CPM+1) |
| H3K27me3 | white → pink | 0–6 log2(CPM+1) |
| Log2FC | blue → white → orange | −2.75–0–2.75 |

---

## Usage

Set working directory to the repository root and run the main script:

```r
# From repository root:
source("single_TE_insertions/scripts/single_TE_insertion_hm.R")
```

Input data are downloaded automatically from the FTP server. Configure analysis parameters at the top of the script under `USER-CONFIGURABLE PARAMETERS`:

```r
SUB_TEs            <- FALSE  # TRUE: plot selected TE family subset; FALSE: all TEs
EUCH_FILT          <- TRUE   # Filter for euchromatic insertions
MIN_NEXT_TE        <- 1500   # Min distance (nt) to neighboring TE
APPLY_H3K4_FILTER  <- TRUE   # Remove TE insertions upon both these conditions
	# 1: under control conditions has a significant amount of H3K4me3 level (>2) 
	# 2: H3K4me3 signal does not change substancially upon piwi knockdown 
FILTER_MODE        <- "union" # "union" or "intersection"
INCLUDE_H3K27      <- FALSE  # Include H3K27me3 panels (moved to supplement)
export_table       <- FALSE  # Export Excel tables of heatmap values
```

---

## Output files

All outputs are written to `output/`:

| File | Description |
|------|-------------|
| `te_grouped_k4_k9_bam.pdf` | Heatmap for BamMut background (GSC-like cells) |
| `te_grouped_k4_k9_mtd.pdf` | Heatmap for MTD background (512c nurse cells) |
| `combined_filter_summary.csv` | Per-insertion H3K4me3 values and filter decisions for both datasets |
| `combined_filtered_elements.txt` | Names of TE insertions removed by H3K4me3 filter |
| `heatmap_creation.log` | Full run log with parameters and summary statistics |
| `*scatter*.pdf` | Scatter plots of H3K9me3 vs H3K4me3 signal (only when `SUB_TEs = TRUE`) |
| `CutnTag_bam_Table_for_heatmap_single_TE_insertions.xlsx` | Matrix values used for heatmap (only when `export_table = TRUE`) |
| `CutnTag_single_TE_insertions_table.xlsx` | Genomic coordinates and metadata of plotted insertions (only when `export_table = TRUE`) |

---

## Scripts

| Script | Description |
|--------|-------------|
| [`scripts/single_TE_insertion_hm.R`](scripts/single_TE_insertion_hm.R) | Main script — data loading, filtering, and heatmap generation |
| [`scripts/cat_hm_helper.R`](scripts/cat_hm_helper.R) | Helper functions: matrix normalization, heatmap building, H3K4me3 filtering logic |
| [`scripts/aesthetics_functions.R`](scripts/aesthetics_functions.R) | Plot style presets (`paper_small`, `paper_large`, `paper_smaller`, `poster`) |

---

## R package dependencies

```r
# Bioconductor
BiocManager::install(c("GenomicRanges", "rtracklayer", "EnrichedHeatmap",
                       "ComplexHeatmap", "GenomeInfoDb"))

# CRAN
install.packages(c("tidyverse", "RColorBrewer", "circlize",
                   "stringr", "ggrastr", "writexl"))
```

---

## Notes

- The `lenfrac` metadata column records the fraction of the TE consensus covered by the nanopore-assembled insertion (reference TEs only). Size information is not available for non-reference insertions detected by ngs_te_mapper2.
- Scatter plots (`SUB_TEs = TRUE`) compare upstream H3K9me3 signal in control against H3K4me3 signal upon *piwi* knockdown, and include log2 fold-change plots for each TE family group.
- H3K27me3 data are excluded from the main figure heatmaps (`INCLUDE_H3K27 = FALSE`) and shown in Supplementary Figure 9.
