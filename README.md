
# spatialCooccur <img src="man/figures/logo.png" align="right" height="138" />

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17498341.svg)](https://doi.org/10.5281/zenodo.17498341)
[![Website](https://img.shields.io/badge/docs-juninamo.github.io%2FspatialCooccur-4a3aa7)](https://juninamo.github.io/spatialCooccur/)

`spatialCooccur` is an R package for analyzing spatial co-occurrence and
neighborhood interactions in spatial transcriptomics data. It quantifies
whether cell types sit together more (or less) often than expected by
chance, maps where they meet, **compares these scores between patient
groups** with statistics that treat the patient as the unit of analysis,
and (experimentally) measures co-localization **directly from transcript
coordinates** without cell segmentation.

**Documentation, tutorials and function reference (searchable):
<https://juninamo.github.io/spatialCooccur/>**

## Installation

``` r
# install.packages("devtools")
devtools::install_github("juninamo/spatialCooccur")
```

## Tutorials

| Tutorial | What it covers |
|---|---|
| [Simulation data: SNA & sCLS](https://juninamo.github.io/spatialCooccur/articles/SNA_tutorial_simulation.html) | Neighborhood enrichment and the local co-localization score on simulated tissue |
| [10x Xenium data: SNA & sCLS](https://juninamo.github.io/spatialCooccur/articles/SNA_tutorial_10Xdata.html) | The same analyses on public Xenium human breast and mouse brain data |
| [Case-control comparison](https://juninamo.github.io/spatialCooccur/articles/case_control_tutorial.html) | Several images per patient: choosing the score, patient-level tests (Wilcoxon / LMM / blocked permutation), pseudoreplication, covariates, power |
| [Segmentation-free co-localization](https://juninamo.github.io/spatialCooccur/articles/segmentation_free_tutorial.html) (experimental) | Cross pair correlation of marker transcripts, a random-feature log-Gaussian Cox process model, case-control testing, Xenium mouse brain |
| [Algorithm reference](https://juninamo.github.io/spatialCooccur/articles/algorithms.html) | The mathematics behind every core function |

The notebooks behind these pages are in [`vignettes/`](vignettes/).

## Features

**Single-sample analysis**

- Simulate spatial layouts with a planted interaction: `generate_sim()`
- Neighborhood enrichment with a label-permutation null:
  `nhood_enrichment()` returns the z-score (evidence) and
  `log2_oe` = log2(observed / expected) (effect size)
- Radius-based co-occurrence ratio: `calc_co_occurrence_for_radius()` /
  `compute_co_occurrence_ratio()`
- Local co-localization: `cooccur_local_oe()` counts A-B pairs around each
  cell, divides by their exact expectation under label permutation and
  smooths with a Gaussian kernel (abundance-adjusted local log2 O/E, optional
  permutation hotspots, O(n k)); the original diffusion sCLS is kept as
  `cooccur_local()`
- Connected interaction spots: `search_interaction_spot()`

**Multi-sample / disease-group comparison**

- Per-image scores for Seurat objects, lists of Seurat objects or plain
  tables: `nhood_enrichment_per_sample()`, `cooccur_ratio_per_sample()`,
  `cooccur_local_per_sample()`, `interaction_spot_per_sample()`;
  `summarize_by_patient()` averages images within patients
- `compare_groups()`: Wilcoxon (exact for small samples), Welch's *t*,
  linear mixed model (`lme4`, Satterthwaite df via `lmerTest`, covariates),
  patient-blocked permutation, and paired designs such as pre- vs
  post-treatment (`method = "signrank"`, within-patient permutation);
  warns on pseudoreplication
- Plots: `plot_group_delta_heatmap()`, `plot_pair_boxplot()`,
  `plot_volcano_groups()`

**Segmentation-free analysis of transcript coordinates (experimental)**

- `read_xenium_transcripts()`, `bin_transcripts()`
- Model-free cross pair correlation of gene sets: `pcf_cross()`,
  `pcf_matrix()`; the relative version is the transcript-level
  counterpart of `log2_oe`
- Random-feature log-Gaussian Cox process factor model after Gundersen,
  Zhang & Engelhardt (AISTATS 2021): `fit_spatial_rff()`,
  `rff_pair_correlation()`
- `colocalization_per_sample()` feeds `compare_groups()`

## How spatialCooccur fits with related tools

Spatial transcriptomics now has excellent tools for discovering structure
directly from the data. spatialCooccur is designed to sit next to them and
focuses on one question: **how strongly do two defined cell populations or
gene programmes co-localize, at which distance, and does this differ between
groups of patients?**

| Tool | Primary focus |
|---|---|
| [FICTURE](https://github.com/seqscope/ficture) (Si *et al.*, *Nat Methods* 2024) | Segmentation-free inference of spatial factors at submicron, pixel-level resolution with a multilayer Dirichlet model; scales to billions of transcripts |
| [punkst](https://github.com/Yichen-Si/punkst) | Scalable toolkit implementing the FICTURE pixel-level factor pipeline and preparing results for visualization |
| [MultiScale_ComplementMacrophage](https://github.com/fanzhanglab/MultiScale_ComplementMacrophage) (Guo *et al.*, in submission) | Spatial neighbourhood-based regression of gene-level associations (Gaussian-kernel neighbourhood exposure, adjusted for self expression and cell density) to define complement-associated niches in RA synovium |
| **spatialCooccur** | Calibrated co-localization between labelled populations (cells, transcripts or factors) and **patient-level comparison across groups** |

What spatialCooccur adds:

- **A calibrated effect size.** Every score is observed / expected under
  label permutation with positions fixed (`log2_oe`, relative pair
  correlation). It is verified on negative controls to stay at zero for any
  number of cell types, cell abundance, image size or cellularity, so values
  can be compared between samples.
- **Distance as an explicit axis.** `pcf_matrix()` returns co-localization as
  a curve over distance for all pairs at once (FFT; 45 pairs over a whole Xenium section in
  about 1.5 minutes), so short-range contact and tissue-scale compartments are
  separated.
- **Where, not only whether.** `cooccur_local_oe()` maps local O/E and
  permutation hotspots for any pair.
- **Case-control and paired designs.** `compare_groups()` treats the patient
  as the unit (exact Wilcoxon, mixed models, blocked permutation, signed-rank
  and sign-flip for pre / post treatment), with power guidance and a
  pseudoreplication warning.
- **Cells, transcripts or factors as input.** The same statistics run on
  segmented cell types, on marker transcripts without segmentation, or on
  any labelled points, for example pixel-level factors from FICTURE / punkst:

  ```r
  # px: pixel-level factor output with coordinates and the top factor
  b <- bin_transcripts(px, bin_size = 4, x_col = "X", y_col = "Y", gene_col = "K1")
  pcf_matrix(b, list(F1 = "1", F2 = "2", F3 = "3"), r_max = 100)
  ```

- **A generative model when needed.** `fit_spatial_rff()` fits a
  random-feature log-Gaussian Cox process (Gundersen, Zhang & Engelhardt,
  AISTATS 2021) that learns spatial length scales and separates cellularity
  from composition.

> **Note for users of versions <= 0.99.1.** Version 0.99.2 fixes the
> permutation null of `nhood_enrichment()` (same-type z-scores were
> inflated) and the diffusion of `cooccur_local()` for `maxnsteps > 1`.
> Results from earlier versions are not directly comparable; see
> [NEWS](NEWS.md).

## Quick start

### 1. Spatial Neighborhood Analysis (SNA)

```r
library(spatialCooccur)
df <- generate_sim(close_ratio = 1, n_types = 15, max_loc = 800, n_cells = 500,
                   test_type = "circle", distance_param = 20, seed = 1234)

res <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = 30,
                        n_perms = 100, seed = 1234, n_jobs = 1)
res$zscore    # evidence: grows with the number of cells
res$log2_oe   # effect size: use this to compare samples
```

### 2. Local co-localization (where do A and B meet?)

```r
rownames(df) <- paste0("cell", seq_len(nrow(df)))
lo <- cooccur_local_oe(df, cluster_x = "cell_type_1", cluster_y = "cell_type_2",
                       radius = 30, n_perms = 99)
attr(lo, "section_log2_oe")   # whole-section effect size
head(lo)                       # per-cell local_log2_oe, p, padj (hotspots)
```

### 3. Comparing patient groups

```r
df_groups <- generate_sim_groups(
  n_samples_per_group = 6, n_images_per_patient = 3,
  group_close_ratio = list(case = 0.8, control = 0.2),
  n_types = 5, n_cells = 400, test_type = "distribute", distance_param = 15
)
per_image <- nhood_enrichment_per_sample(
  df_groups, sample_key = "sample_id", group_key = "group",
  cluster_key = "cell_type", patient_key = "patient",
  neighbors.k = 20, n_perms = 100, n_jobs = 1
)

# images are nested in patients: use a mixed model (or aggregate first
# with summarize_by_patient() and use method = "wilcox")
res <- compare_groups(per_image, value = "log2_oe", method = "lmm",
                      patient_key = "patient", ref_group = "control",
                      symmetric = TRUE)
head(res)
plot_volcano_groups(res, label_top = 5)
```

### 4. Segmentation-free co-localization (experimental)

```r
tx <- read_xenium_transcripts("path/to/xenium_outs", genes = unlist(marker_sets))
b  <- bin_transcripts(tx, bin_size = 4)
pcf_matrix(b, marker_sets, r_max = 100)   # log_g and log_g_rel for every pair
```

## Figures from the paper

<!-- 
## Citation 
Jun Inamo, Roselyn Fierkens, Michael R CLay, Anna Helena Jonsson, Clara Lin, Kari Hayes, Nathan Rogers, Heather Leach, Kentaro Yomogida. Subcellular spatial transcriptomics reveals immune–stromal crosstalk within the synovium of patients with juvenile idiopathic arthritis. [*bioRxiv*](https://www.biorxiv.org/XX), doi:[https://doi.org/XX](https://doi.org/XX)

- Study design and identified cell clusters in JIA-synovium

<kbd>
<img src="man/figures/Figure1.png" width="800" align="center">
</kbd>

&nbsp;&nbsp;
-->

- Spatial Neighborhood Analysis

<kbd>
<img src="man/figures/Figure2.png" width="800" align="center">
</kbd>

&nbsp;&nbsp;

- Spatial Co-localization Score

<kbd>
<img src="man/figures/Figure3.png" width="800" align="center">
</kbd>

&nbsp;&nbsp;

## 📝 Citation 
Jun Inamo, Roselyn Fierkens, Michael R. Clay, Anna Helena Jonsson, Clara Lin, Kari Hayes, Nathan Rogers, Heather Leach, Kentaro Yomogida. Spatial transcriptomics reveals immune–stromal crosstalk within the synovium of patients with juvenile idiopathic arthritis. [*JCI Insight* 2026;11(1):e198074](https://doi.org/10.1172/jci.insight.198074). doi:[10.1172/jci.insight.198074](https://doi.org/10.1172/jci.insight.198074)

## Contact
For questions or issues, please open a GitHub issue or contact:

**Name:** Jun Inamo  
**Email:** juninamo@keio.jp  
**Affiliation:** Department of Microbiology and Immunology, Keio University School of Medicine

<!-- 
## Acknowledgments
This work was supported by the Uehara Memorial Foundation Postdoctoral Fellowship, a Grant-in-Aid for Japan Society for the Promotion of Science Overseas Research Fellows, the Mochida Memorial Foundation for Medical and Pharmaceutical Research (to J.I.), K08DK128544 (K.Y.). 
-->
&nbsp;&nbsp;

## License
This repository is provided under the MIT License.

