# spatialCooccur

[![DOI](https://zenodo.org/badge/960705353.svg)](https://doi.org/10.5281/zenodo.17498341)

`spatialCooccur` is an R package for analyzing spatial co-occurrence and
neighborhood interactions in spatial transcriptomics data. It is built
around Seurat objects and provides tools to compute co-occurrence
enrichment, perform permutation-based tests, visualize local interaction
scores, and **compare scores between disease groups** across multiple
samples.

## Installation

You can install the development version from GitHub using:

``` r

# install.packages("devtools")
devtools::install_github("juninamo/spatialCooccur")
```

## Features

**Single-sample analysis**

- Simulate spatial transcriptomic layouts with
  [`generate_sim()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim.md)
- Calculate neighborhood co-occurrence enrichment with
  [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
  (permutation-based z-score)
- Compute radius-based co-occurrence ratio with
  [`calc_co_occurrence_for_radius()`](https://juninamo.github.io/spatialCooccur/reference/calc_co_occurrence_for_radius.md)
  /
  [`compute_co_occurrence_ratio()`](https://juninamo.github.io/spatialCooccur/reference/compute_co_occurrence_ratio.md)
- Identify local interaction zones using
  [`cooccur_local()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.md)
- Detect connected interaction spots with
  [`search_interaction_spot()`](https://juninamo.github.io/spatialCooccur/reference/search_interaction_spot.md)

**Multi-sample / disease-group comparison**

- Generate multi-sample group-structured simulations with
  [`generate_sim_groups()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim_groups.md)
- Compute per-sample scores with
  [`nhood_enrichment_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment_per_sample.md),
  [`cooccur_ratio_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_ratio_per_sample.md),
  [`cooccur_local_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local_per_sample.md),
  or
  [`interaction_spot_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/interaction_spot_per_sample.md)
  (supports Seurat / list of Seurat / data.frame input, image- or
  patient-level aggregation)
- Test cluster pairs between groups with
  [`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
  — Wilcoxon / Welch’s *t* / linear mixed model (`lme4`) /
  patient-blocked permutation
- Visualize with
  [`plot_group_delta_heatmap()`](https://juninamo.github.io/spatialCooccur/reference/plot_group_delta_heatmap.md),
  [`plot_pair_boxplot()`](https://juninamo.github.io/spatialCooccur/reference/plot_pair_boxplot.md),
  [`plot_volcano_groups()`](https://juninamo.github.io/spatialCooccur/reference/plot_volcano_groups.md)

**Works with Seurat spatial objects out of the box.**

**Vignettes**

- `vignette("disease_comparison", "spatialCooccur")` — end-to-end
  disease-group comparison workflow with worked example
- [`vignette("algorithms", "spatialCooccur")`](https://juninamo.github.io/spatialCooccur/articles/algorithms.md)
  — mathematical reference for every core function
- `vignettes/case_control_tutorial.ipynb` — case-control study with
  several images per patient: choosing the score, patient-level testing
  (Wilcoxon / LMM / blocked permutation), pseudoreplication, covariates,
  and power

### 1. Spatial Neighborhood Analysis (SNA)

To simulate spatial transcriptomic data and perform neighborhood
enrichment analysis:

``` r

df = generate_sim(close_ratio = 1, n_types = 15, max_loc = 800, n_cells = 500, test_type = "circle", distance_param = 20, seed=1234)

# Run neighborhood enrichment analysis
nhood_enrichment_res <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = 30, n_perms = 100, seed = 1234, n_jobs = 4)
nhood_enrichment_res$zscore
```

### 2. Spatial Co-localization Score (sCLS)

To compute co-localization scores for cell interactions:

``` r

cooccur_local_df <- cooccur_local(df, cluster_x = "cell_type_1", cluster_y = "cell_type_2", neighbors.k = 30, radius = 30)
summary(cooccur_local_df)
```

### 3. Disease-group comparison

To compare a spatial co-occurrence score between disease groups across
multiple samples, compute per-sample scores and then test each cluster
pair. The same `*_per_sample()` +
[`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
pattern works for neighborhood enrichment z-score, radius-based ratio,
local co-occurrence score, and interaction-spot counts.

``` r

# Simulate two groups, 3 samples each
df_groups <- generate_sim_groups(
  n_samples_per_group = 3,
  group_close_ratio = list(disease = 0.8, control = 0.2),
  n_types = 5, n_cells = 400, test_type = "distribute",
  distance_param = 15
)

# Per-sample z-scores (one row per sample x cluster_i x cluster_j)
per_sample <- nhood_enrichment_per_sample(
  df_groups, sample_key = "sample_id", group_key = "group",
  cluster_key = "cell_type", patient_key = "patient",
  neighbors.k = 20, n_perms = 100
)

# Group comparison: method = "wilcox" (default) | "t" | "lmm" | "perm"
# - "lmm"  uses lme4::lmer(value ~ group + (1 | patient))
# - "perm" runs a group-label permutation test, blocked by patient
#          when patient_key is supplied
# Compare log2(observed / expected) rather than the z-score: z-scores grow
# with the number of cells per image, log2_oe does not.
res <- compare_groups(
  per_sample, value = "log2_oe",
  method = "wilcox", ref_group = "control", symmetric = TRUE
)
head(res)

# Visualize: per-pair effect heatmap, volcano, and per-sample boxplot
plot_group_delta_heatmap(res)
plot_volcano_groups(res, label_top = 5)
plot_pair_boxplot(
  per_sample, value = "log2_oe",
  pairs = data.frame(cluster_i = "cell_type_1",
                     cluster_j = "cell_type_2"),
  add_p = TRUE, ref_group = "control"
)
```

See `vignette("disease_comparison", "spatialCooccur")` for the full
workflow, sanity-check heatmaps, and the algorithm description (math)
for each
[`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
method.

- Spatial Neighborhood Analysis

![](reference/figures/Figure2.png)

  

- Spatial Co-localization Score

![](reference/figures/Figure3.png)

  

## 📝 Citation

Jun Inamo, Roselyn Fierkens, Michael R. Clay, Anna Helena Jonsson, Clara
Lin, Kari Hayes, Nathan Rogers, Heather Leach, Kentaro Yomogida. Spatial
transcriptomics reveals immune–stromal crosstalk within the synovium of
patients with juvenile idiopathic arthritis. [*JCI Insight*
2026;11(1):e198074](https://doi.org/10.1172/jci.insight.198074).
<doi:%5B10.1172/jci.insight.198074>\](<https://doi.org/10.1172/jci.insight.198074>)

## Contact

For questions or issues related to this tutorial, please contact;

**Name:** Jun Inamo  
**Email:** <juninamo@keio.jp>  
**Affiliation:** Department of Microbiology and Immunology, Keio
University School of Medicine

The data presented in the paper (spatial transcriptome data from
JIA-synovoum) was generated by the [Yomogida
lab](https://www.yomogidalab.com/).

  

## License

This repository is provided under the MIT License.
