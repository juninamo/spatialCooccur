# cohalu 0.99.3

## Package renamed

* The package is renamed from `spatialCooccur` to `cohalu` (COHALU:
  CO-localization, Hotspots And sample-Level Units; read "koharu"). Use
  `library(cohalu)`; function names and arguments are unchanged. The internal
  class `spatialCooccurSample` is now `cohaluSample`.

## Changes that affect results

* `nhood_enrichment()`, `nhood_enrichment.Seurat()` and
  `nhood_enrichment_per_sample()`: the default is now `transformation = FALSE`
  (every kNN link counts once, as in squidpy). The previous default weighted
  each link from cell u by 1 / (1 + d_u), d_u = number of cells that chose u.
  With a kNN graph every cell sends k links at any density, so this was not a
  density correction; in simulations it kept the calibration (false positives
  4-5% per pair, 2-7% family-wise) but lowered power (planted pair, 1,500
  cells: 14% with the weighting vs 28% without). Use `transformation = TRUE`
  to reproduce earlier results. The docs now describe what the weighting does.
* `nhood_enrichment()`: `log2_oe` is now centred on the label shuffles
  (the mean of the same log ratio over the shuffles is subtracted), so it is
  0 on average without interaction for any number of cells. The log of a
  ratio of small counts was biased below 0 for rare cell types (about -0.03
  for 25 cell types of 60 cells, -0.05 for pairs with < 50 cells). The
  previous value is returned as `log2_oe_raw`.
* `generate_sim(test_type = "distribute")` no longer places relocated cells
  outside the tissue. Such cells sat in empty space where their k nearest
  neighbours reached their partner cell even at 100 um, so co-localization
  leaked to long planted distances.
* `compare_groups()`: new `unit = c("patient", "image")`, default
  `"patient"`. With `patient_key` and `method = "wilcox"` or `"t"`, images
  are averaged within patient before testing, so the patient is the unit of
  analysis. Previously these tests used image-level rows (pseudoreplication,
  13.8% false positives at 5% in simulations) unless the data were first
  passed through `summarize_by_patient()`; use `unit = "image"` for the old
  behaviour.
* `compare_groups()` now defaults to `value = "log2_oe"` and uses a
  `patient` column automatically when `patient_key` is not given and
  patients have several images.
* `compare_groups(symmetric = TRUE)` averages the (i, j) and (j, i) values
  of each sample instead of keeping only the `cluster_i <= cluster_j` row
  (degree-normalised neighbourhood scores are slightly directional).
* `summarize_by_patient()` keeps distances (`r`) of
  `colocalization_per_sample()` output separate instead of averaging them;
  new `pair_keys` argument.

## New features

* Residual random-feature model (experimental): `fit_spatial_rff()` accepts
  a per-bin, per-gene log `offset` matrix describing structure that is
  already known (cell-type composition, domains, or an embedding such as
  PCA, Harmony or SCIGMA; built with the new `rff_offset()`), so that the
  factors capture only the remaining spatially coherent variation. New
  `ard` argument: a group penalty on each factor's loadings that shrinks
  unneeded factors, reported as `factor_strength`. New `rff_factor_test()`:
  parametric-bootstrap p-value per factor against the largest factor fitted
  to null data (family-wise over factors). The test statistic is the program
  strength (loading norm after removing the loading shared by all genes):
  factors that move all genes together are cellularity, not gene programs,
  and the raw loading norm called such factors significant in simulations
  without any program. `fit_spatial_rff()` also returns `program_strength`.
* `rff_fields()` (experimental): evaluates the latent fields of a
  `fit_spatial_rff()` fit (cellularity and the `K` factors, on the
  unit-variance scale of `fit$field_grid`) at any coordinates - cell
  centroids, transcripts or other points - and averages them per cell with
  `by = "cell_id"`. `fit_spatial_rff()` now also stores the coordinate
  centre used for the random features (`center`); older fits still work.
* `associate_continuous()`: association of per-image / per-patient
  co-localization with a continuous clinical variable (CRP, disease
  activity, age): Spearman on patient means (default), linear model with
  covariates, mixed model on images, or permutation; BH over pairs.
* Unsupervised transcript-level co-localization (experimental):
  `colocalization_gene_matrix()` (gene x gene log2 O/E of transcript pairs
  within a radius, label-shuffling expectation in closed form, one FFT per
  gene), `colocalization_modules()` (clusters co-localizing genes),
  `module_enrichment()` (hypergeometric test of any gene sets, e.g. pathways
  or cell-type markers) and `module_enrichr()` (enrichR wrapper).
* `nhood_enrichment()` returns a within-sample test per unordered pair:
  `pvalue` (normal, from the shuffles), `padj` (Westfall-Young max-T,
  family-wise error rate, calibrated for any number of cell types) and
  `padj_bh`. New `plot_nhood_heatmap()` draws `log2_oe` with significance
  stars.
* `nhood_enrichment()` also returns directional statistics (row = centre
  cell type, column = neighbour type; not symmetric), because the pair-level
  `log2_oe` is symmetric by construction and cannot tell "A is surrounded by
  B" from "B is surrounded by A": `contact` (share of centre cells with at
  least one neighbour of the other type) and `dominance` (share of centre
  cells whose neighbours are at least half of the other type), each with
  `*_expected`, centred `*_log2_oe`, `*_pvalue` and max-T `*_padj` over all
  ordered pairs. `plot_nhood_heatmap(value = "dominance_log2_oe")` draws them
  without symmetrising; the pair-level heatmap is labelled as the average of
  both directions.
* `plot_nhood_heatmap()`: new default `triangle = "auto"` draws symmetric
  pair-level values once (lower triangle with the diagonal) and directional
  values in full; `triangle = "full"` restores the previous layout.
* Tutorials use `log2_oe` with the max-T adjusted `padj`; new sections on
  `cooccur_local_oe()`, `associate_continuous()` and gene-level modules. The
  algorithm reference covers all current methods.

# cohalu 0.99.2

## Bug fixes that change results

* `nhood_enrichment()` / `permute_clusters()`: the permutation null now applies
  a single label permutation to rows and columns. Previously rows and columns
  were shuffled independently, which under spatial randomness inflated
  same-type z-scores (about +11) and biased different-type z-scores (about -2).
  z-scores are not comparable with versions <= 0.99.1.
* `cooccur_local()` / `cooccur_local.Seurat()`: diffusion now iterates from the
  previous step, and the kurtosis early-stop rule compares consecutive steps.
  Previously every step restarted from the raw indicator, so any
  `maxnsteps >= 1` gave the one-step result. Results with `maxnsteps <= 1`
  are unchanged.

## Other bug fixes

* `nhood_enrichment()` with `n_jobs > 1`: the workers did not load the Matrix
  methods, so every worker failed and the permutations silently fell back to
  sequential (no speed-up). The workers now load Matrix, and a fallback to
  sequential is reported as a warning.
* `nhood_enrichment()` is now reproducible for a given `seed`, both
  sequentially and with `n_jobs > 1` (worker RNG streams are seeded).
* `compare_groups()` and `nhood_enrichment()` no longer overwrite the caller's
  RNG state.
* Seurat input to the `*_per_sample()` helpers matches images to `meta.data`
  by cell name, so `sample_key` values need not equal the image names.
* `interaction_spot_per_sample()` returns `NA` (not 0 spots) when the spot
  search fails, and no longer requires a `cell` column in `meta.data`.
* Absent cell types give `NA` scores instead of silently dropped rows.

## New features

* `cooccur_local_oe()`: abundance-adjusted local co-localization. For every
  cell, the number of cluster_x-cluster_y pairs within `radius` is divided by
  its exact expectation under label permutation (closed form), smoothed with
  a Gaussian kernel of explicit width, with optional permutation hotspot
  p-values (O(n k) per permutation). `cooccur_local_per_sample()` gains the
  corresponding `log2_oe` summary, recommended for group comparison: the
  mean diffusion sCLS is unchanged by diffusion (mass-conserving) and grows
  with cell-type abundance.

* **Experimental segmentation-free analysis** of transcript coordinates
  (e.g. Xenium `transcripts.parquet`): `bin_transcripts()`, model-free cross
  pair correlation of gene sets `pcf_cross()` (the relative version is the
  label-permutation O/E, the continuous analogue of `log2_oe`), a
  random-feature log-Gaussian Cox process factor model `fit_spatial_rff()`
  after Gundersen, Zhang & Engelhardt (AISTATS 2021) with model-based
  `rff_pair_correlation()`, and `colocalization_per_sample()` whose output
  goes into `compare_groups()` with `pair_keys = c("cluster_i", "cluster_j", "r")`.
  Simulators `simulate_transcripts()` / `simulate_transcripts_groups()` and
  `lgcp_true_pair_correlation()` support validation.
  `read_xenium_transcripts()` reads Xenium transcript tables (binary gene
  names in older outputs, qv filter, gene filtering inside Arrow for 5K
  panels) and `pcf_matrix()` computes all gene-set pairs with cached FFTs.
  Any labelled point set can be analysed the same way, e.g. pixel-level
  factors from FICTURE / punkst (`bin_transcripts(gene_col = "K1")`).
* `compare_groups()` supports paired / repeated-measures designs
  (e.g. pre- vs post-treatment): new `method = "signrank"`, and
  `method = "perm"` now permutes labels within patients when patients appear
  in both groups (the previous between-patient permutation was invalid for
  such designs).

* `nhood_enrichment()` also returns `expected` (permutation mean) and
  `log2_oe` (log2 observed / expected), an effect size that does not grow with
  the number of cells. Recommended for between-group comparison.
* `*_per_sample()` outputs gain `n_cells`, `n_i`, `n_j` (abundance of the two
  cell types).
* `compare_groups()`: exact Wilcoxon p-values for small samples; LMM p-values
  with Satterthwaite df (via lmerTest) instead of Wald z; `covariates`;
  `symmetric`; `min_n_per_group`; exact enumeration for small permutation
  tests; a warning on pseudoreplication.
* `summarize_by_patient()` aggregates image-level results to patients.
* `generate_sim_groups()` supports several images per patient
  (`n_images_per_patient`, `within_patient_noise`) and group-specific
  `n_cells`.
* New tutorial: case-control comparison (`vignettes/case_control_tutorial.ipynb`).
* Faster counting with sparse matrix products.
