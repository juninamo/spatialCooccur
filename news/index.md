# Changelog

## cohalu 0.99.3

### Package renamed

- The package is renamed from `spatialCooccur` to `cohalu` (COHALU:
  CO-localization, Hotspots And sample-Level Units; read “koharu”). Use
  [`library(cohalu)`](https://juninamo.github.io/cohalu/); function
  names and arguments are unchanged. The internal class
  `spatialCooccurSample` is now `cohaluSample`.

### Changes that affect results

- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md):
  `log2_oe` is now centred on the label shuffles (the mean of the same
  log ratio over the shuffles is subtracted), so it is 0 on average
  without interaction for any number of cells. The log of a ratio of
  small counts was biased below 0 for rare cell types (about -0.03 for
  25 cell types of 60 cells, -0.05 for pairs with \< 50 cells). The
  previous value is returned as `log2_oe_raw`.
- `generate_sim(test_type = "distribute")` no longer places relocated
  cells outside the tissue. Such cells sat in empty space where their k
  nearest neighbours reached their partner cell even at 100 um, so
  co-localization leaked to long planted distances.
- [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md):
  new `unit = c("patient", "image")`, default `"patient"`. With
  `patient_key` and `method = "wilcox"` or `"t"`, images are averaged
  within patient before testing, so the patient is the unit of analysis.
  Previously these tests used image-level rows (pseudoreplication, 13.8%
  false positives at 5% in simulations) unless the data were first
  passed through
  [`summarize_by_patient()`](https://juninamo.github.io/cohalu/reference/summarize_by_patient.md);
  use `unit = "image"` for the old behaviour.
- [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md)
  now defaults to `value = "log2_oe"` and uses a `patient` column
  automatically when `patient_key` is not given and patients have
  several images.
- `compare_groups(symmetric = TRUE)` averages the (i, j) and (j, i)
  values of each sample instead of keeping only the
  `cluster_i <= cluster_j` row (degree-normalised neighbourhood scores
  are slightly directional).
- [`summarize_by_patient()`](https://juninamo.github.io/cohalu/reference/summarize_by_patient.md)
  keeps distances (`r`) of
  [`colocalization_per_sample()`](https://juninamo.github.io/cohalu/reference/colocalization_per_sample.md)
  output separate instead of averaging them; new `pair_keys` argument.

### New features

- Residual random-feature model (experimental):
  [`fit_spatial_rff()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff.md)
  accepts a per-bin, per-gene log `offset` matrix describing structure
  that is already known (cell-type composition, domains, or an embedding
  such as PCA, Harmony or SCIGMA; built with the new
  [`rff_offset()`](https://juninamo.github.io/cohalu/reference/rff_offset.md)),
  so that the factors capture only the remaining spatially coherent
  variation. New `ard` argument: a group penalty on each factor’s
  loadings that shrinks unneeded factors, reported as `factor_strength`.
  New
  [`rff_factor_test()`](https://juninamo.github.io/cohalu/reference/rff_factor_test.md):
  parametric-bootstrap p-value per factor against the largest factor
  fitted to null data (family-wise over factors).
- [`rff_fields()`](https://juninamo.github.io/cohalu/reference/rff_fields.md)
  (experimental): evaluates the latent fields of a
  [`fit_spatial_rff()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff.md)
  fit (cellularity and the `K` factors, on the unit-variance scale of
  `fit$field_grid`) at any coordinates - cell centroids, transcripts or
  other points - and averages them per cell with `by = "cell_id"`.
  [`fit_spatial_rff()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff.md)
  now also stores the coordinate centre used for the random features
  (`center`); older fits still work.
- [`associate_continuous()`](https://juninamo.github.io/cohalu/reference/associate_continuous.md):
  association of per-image / per-patient co-localization with a
  continuous clinical variable (CRP, disease activity, age): Spearman on
  patient means (default), linear model with covariates, mixed model on
  images, or permutation; BH over pairs.
- Unsupervised transcript-level co-localization (experimental):
  [`colocalization_gene_matrix()`](https://juninamo.github.io/cohalu/reference/colocalization_gene_matrix.md)
  (gene x gene log2 O/E of transcript pairs within a radius,
  label-shuffling expectation in closed form, one FFT per gene),
  [`colocalization_modules()`](https://juninamo.github.io/cohalu/reference/colocalization_modules.md)
  (clusters co-localizing genes),
  [`module_enrichment()`](https://juninamo.github.io/cohalu/reference/module_enrichment.md)
  (hypergeometric test of any gene sets, e.g. pathways or cell-type
  markers) and
  [`module_enrichr()`](https://juninamo.github.io/cohalu/reference/module_enrichr.md)
  (enrichR wrapper).
- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  returns a within-sample test per unordered pair: `pvalue` (normal,
  from the shuffles), `padj` (Westfall-Young max-T, family-wise error
  rate, calibrated for any number of cell types) and `padj_bh`. New
  [`plot_nhood_heatmap()`](https://juninamo.github.io/cohalu/reference/plot_nhood_heatmap.md)
  draws `log2_oe` with significance stars.
- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  also returns directional statistics (row = centre cell type, column =
  neighbour type; not symmetric), because the pair-level `log2_oe` is
  symmetric by construction and cannot tell “A is surrounded by B” from
  “B is surrounded by A”: `contact` (share of centre cells with at least
  one neighbour of the other type) and `dominance` (share of centre
  cells whose neighbours are at least half of the other type), each with
  `*_expected`, centred `*_log2_oe`, `*_pvalue` and max-T `*_padj` over
  all ordered pairs. `plot_nhood_heatmap(value = "dominance_log2_oe")`
  draws them without symmetrising; the pair-level heatmap is labelled as
  the average of both directions.
- [`plot_nhood_heatmap()`](https://juninamo.github.io/cohalu/reference/plot_nhood_heatmap.md):
  new default `triangle = "auto"` draws symmetric pair-level values once
  (lower triangle with the diagonal) and directional values in full;
  `triangle = "full"` restores the previous layout.
- Tutorials use `log2_oe` with the max-T adjusted `padj`; new sections
  on
  [`cooccur_local_oe()`](https://juninamo.github.io/cohalu/reference/cooccur_local_oe.md),
  [`associate_continuous()`](https://juninamo.github.io/cohalu/reference/associate_continuous.md)
  and gene-level modules. The algorithm reference covers all current
  methods.

## cohalu 0.99.2

### Bug fixes that change results

- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  /
  [`permute_clusters()`](https://juninamo.github.io/cohalu/reference/permute_clusters.md):
  the permutation null now applies a single label permutation to rows
  and columns. Previously rows and columns were shuffled independently,
  which under spatial randomness inflated same-type z-scores (about +11)
  and biased different-type z-scores (about -2). z-scores are not
  comparable with versions \<= 0.99.1.
- [`cooccur_local()`](https://juninamo.github.io/cohalu/reference/cooccur_local.md)
  /
  [`cooccur_local.Seurat()`](https://juninamo.github.io/cohalu/reference/cooccur_local.Seurat.md):
  diffusion now iterates from the previous step, and the kurtosis
  early-stop rule compares consecutive steps. Previously every step
  restarted from the raw indicator, so any `maxnsteps >= 1` gave the
  one-step result. Results with `maxnsteps <= 1` are unchanged.

### Other bug fixes

- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  is now reproducible for a given `seed`, both sequentially and with
  `n_jobs > 1` (worker RNG streams are seeded).
- [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md)
  and
  [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  no longer overwrite the caller’s RNG state.
- Seurat input to the `*_per_sample()` helpers matches images to
  `meta.data` by cell name, so `sample_key` values need not equal the
  image names.
- [`interaction_spot_per_sample()`](https://juninamo.github.io/cohalu/reference/interaction_spot_per_sample.md)
  returns `NA` (not 0 spots) when the spot search fails, and no longer
  requires a `cell` column in `meta.data`.
- Absent cell types give `NA` scores instead of silently dropped rows.

### New features

- [`cooccur_local_oe()`](https://juninamo.github.io/cohalu/reference/cooccur_local_oe.md):
  abundance-adjusted local co-localization. For every cell, the number
  of cluster_x-cluster_y pairs within `radius` is divided by its exact
  expectation under label permutation (closed form), smoothed with a
  Gaussian kernel of explicit width, with optional permutation hotspot
  p-values (O(n k) per permutation).
  [`cooccur_local_per_sample()`](https://juninamo.github.io/cohalu/reference/cooccur_local_per_sample.md)
  gains the corresponding `log2_oe` summary, recommended for group
  comparison: the mean diffusion sCLS is unchanged by diffusion
  (mass-conserving) and grows with cell-type abundance.

- **Experimental segmentation-free analysis** of transcript coordinates
  (e.g. Xenium `transcripts.parquet`):
  [`bin_transcripts()`](https://juninamo.github.io/cohalu/reference/bin_transcripts.md),
  model-free cross pair correlation of gene sets
  [`pcf_cross()`](https://juninamo.github.io/cohalu/reference/pcf_cross.md)
  (the relative version is the label-permutation O/E, the continuous
  analogue of `log2_oe`), a random-feature log-Gaussian Cox process
  factor model
  [`fit_spatial_rff()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff.md)
  after Gundersen, Zhang & Engelhardt (AISTATS 2021) with model-based
  [`rff_pair_correlation()`](https://juninamo.github.io/cohalu/reference/rff_pair_correlation.md),
  and
  [`colocalization_per_sample()`](https://juninamo.github.io/cohalu/reference/colocalization_per_sample.md)
  whose output goes into
  [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md)
  with `pair_keys = c("cluster_i", "cluster_j", "r")`. Simulators
  [`simulate_transcripts()`](https://juninamo.github.io/cohalu/reference/simulate_transcripts.md)
  /
  [`simulate_transcripts_groups()`](https://juninamo.github.io/cohalu/reference/simulate_transcripts_groups.md)
  and
  [`lgcp_true_pair_correlation()`](https://juninamo.github.io/cohalu/reference/lgcp_true_pair_correlation.md)
  support validation.
  [`read_xenium_transcripts()`](https://juninamo.github.io/cohalu/reference/read_xenium_transcripts.md)
  reads Xenium transcript tables (binary gene names in older outputs, qv
  filter, gene filtering inside Arrow for 5K panels) and
  [`pcf_matrix()`](https://juninamo.github.io/cohalu/reference/pcf_matrix.md)
  computes all gene-set pairs with cached FFTs. Any labelled point set
  can be analysed the same way, e.g. pixel-level factors from FICTURE /
  punkst (`bin_transcripts(gene_col = "K1")`).

- [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md)
  supports paired / repeated-measures designs (e.g. pre- vs
  post-treatment): new `method = "signrank"`, and `method = "perm"` now
  permutes labels within patients when patients appear in both groups
  (the previous between-patient permutation was invalid for such
  designs).

- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  also returns `expected` (permutation mean) and `log2_oe` (log2
  observed / expected), an effect size that does not grow with the
  number of cells. Recommended for between-group comparison.

- `*_per_sample()` outputs gain `n_cells`, `n_i`, `n_j` (abundance of
  the two cell types).

- [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md):
  exact Wilcoxon p-values for small samples; LMM p-values with
  Satterthwaite df (via lmerTest) instead of Wald z; `covariates`;
  `symmetric`; `min_n_per_group`; exact enumeration for small
  permutation tests; a warning on pseudoreplication.

- [`summarize_by_patient()`](https://juninamo.github.io/cohalu/reference/summarize_by_patient.md)
  aggregates image-level results to patients.

- [`generate_sim_groups()`](https://juninamo.github.io/cohalu/reference/generate_sim_groups.md)
  supports several images per patient (`n_images_per_patient`,
  `within_patient_noise`) and group-specific `n_cells`.

- New tutorial: case-control comparison
  (`vignettes/case_control_tutorial.ipynb`).

- Faster counting with sparse matrix products.
