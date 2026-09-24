# Package index

## Simulation

Generate synthetic tissues with a planted cell-cell interaction, and
convert them to Seurat objects.

- [`generate_sim()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim.md)
  : Simulate Spatial Coordinates and Cell Types
- [`generate_sim_groups()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim_groups.md)
  : Generate multi-sample simulated data with disease-group structure
- [`sim_to_seurat()`](https://juninamo.github.io/spatialCooccur/reference/sim_to_seurat.md)
  : Convert simulated cells to a Seurat object with one FOV per sample

## Single-sample analysis

Spatial neighborhood analysis (SNA), spatial co-localization score
(sCLS), radius-based co-occurrence and interaction spots for one tissue.

- [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
  : Neighborhood Enrichment (Generic method)
- [`nhood_enrichment.Seurat()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.Seurat.md)
  : Neighborhood Enrichment (Seurat Method)
- [`cooccur_local()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.md)
  : Local Co-occurrence Score (Generic method)
- [`cooccur_local.Seurat()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.Seurat.md)
  : Local Co-occurrence Score (Seurat Method)
- [`calc_co_occurrence_for_radius()`](https://juninamo.github.io/spatialCooccur/reference/calc_co_occurrence_for_radius.md)
  : Calculate Co-occurrence Matrix for a Given Radius
- [`compute_co_occurrence_ratio()`](https://juninamo.github.io/spatialCooccur/reference/compute_co_occurrence_ratio.md)
  : Compute Enrichment Ratios from Count Matrix
- [`search_interaction_spot()`](https://juninamo.github.io/spatialCooccur/reference/search_interaction_spot.md)
  : Search for Spatial Interaction Spots

## Multi-sample scores

Compute a score per image (or patient) and return a tidy table ready for
group comparison.

- [`build_sample_design()`](https://juninamo.github.io/spatialCooccur/reference/build_sample_design.md)
  : Build a sample design table for disease-group comparisons
- [`nhood_enrichment_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment_per_sample.md)
  : Per-sample neighborhood enrichment
- [`cooccur_ratio_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_ratio_per_sample.md)
  : Per-sample radius-based co-occurrence ratio
- [`cooccur_local_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local_per_sample.md)
  : Per-sample local co-occurrence summary
- [`interaction_spot_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/interaction_spot_per_sample.md)
  : Per-sample interaction-spot summary
- [`summarize_by_patient()`](https://juninamo.github.io/spatialCooccur/reference/summarize_by_patient.md)
  : Aggregate per-image scores to one row per patient

## Segmentation-free analysis (experimental)

Work directly on transcript coordinates: bin transcripts, estimate the
cross pair correlation of gene sets model-free, or fit a random-feature
log-Gaussian Cox process factor model (after Gundersen et al. 2021).

- [`read_xenium_transcripts()`](https://juninamo.github.io/spatialCooccur/reference/read_xenium_transcripts.md)
  : Read a Xenium transcript table
- [`bin_transcripts()`](https://juninamo.github.io/spatialCooccur/reference/bin_transcripts.md)
  : Bin transcript coordinates on a square grid
- [`pcf_cross()`](https://juninamo.github.io/spatialCooccur/reference/pcf_cross.md)
  : Empirical cross pair correlation of two gene sets
  (segmentation-free)
- [`fit_spatial_rff()`](https://juninamo.github.io/spatialCooccur/reference/fit_spatial_rff.md)
  : Fit a random-feature spatial factor model to binned transcripts
- [`rff_pair_correlation()`](https://juninamo.github.io/spatialCooccur/reference/rff_pair_correlation.md)
  : Model-based cross pair correlation from a fitted random-feature LGCP
- [`colocalization_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/colocalization_per_sample.md)
  : Per-sample segmentation-free co-localization of gene sets
- [`simulate_transcripts()`](https://juninamo.github.io/spatialCooccur/reference/simulate_transcripts.md)
  : Simulate transcript coordinates from a multivariate log-Gaussian Cox
  process
- [`simulate_transcripts_groups()`](https://juninamo.github.io/spatialCooccur/reference/simulate_transcripts_groups.md)
  : Simulate a multi-sample case-control transcript dataset
- [`lgcp_true_pair_correlation()`](https://juninamo.github.io/spatialCooccur/reference/lgcp_true_pair_correlation.md)
  : True cross pair correlation of a simulated LGCP

## Group comparison

Test every cell-type pair between two groups.

- [`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
  : Compare disease groups across samples

## Visualization

- [`plot_group_delta_heatmap()`](https://juninamo.github.io/spatialCooccur/reference/plot_group_delta_heatmap.md)
  : Heatmap of the disease-group effect for every cluster pair
- [`plot_pair_boxplot()`](https://juninamo.github.io/spatialCooccur/reference/plot_pair_boxplot.md)
  : Per-sample boxplot of one (or several) cluster pair(s) across groups
- [`plot_volcano_groups()`](https://juninamo.github.io/spatialCooccur/reference/plot_volcano_groups.md)
  : Volcano plot of a compare_groups result
- [`manual_colors`](https://juninamo.github.io/spatialCooccur/reference/manual_colors.md)
  : Default Manual Colors for Clusters

## Low-level helpers

- [`compute_count()`](https://juninamo.github.io/spatialCooccur/reference/compute_count.md)
  : Compute Co-occurrence Count Matrix
- [`permute_clusters()`](https://juninamo.github.io/spatialCooccur/reference/permute_clusters.md)
  : Permute Cluster Assignments and Recompute Counts
