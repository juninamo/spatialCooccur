# Changelog

## spatialCooccur 0.99.2

### Bug fixes that change results

- [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
  /
  [`permute_clusters()`](https://juninamo.github.io/spatialCooccur/reference/permute_clusters.md):
  the permutation null now applies a single label permutation to rows
  and columns. Previously rows and columns were shuffled independently,
  which under spatial randomness inflated same-type z-scores (about +11)
  and biased different-type z-scores (about -2). z-scores are not
  comparable with versions \<= 0.99.1.
- [`cooccur_local()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.md)
  /
  [`cooccur_local.Seurat()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.Seurat.md):
  diffusion now iterates from the previous step, and the kurtosis
  early-stop rule compares consecutive steps. Previously every step
  restarted from the raw indicator, so any `maxnsteps >= 1` gave the
  one-step result. Results with `maxnsteps <= 1` are unchanged.

### Other bug fixes

- [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
  is now reproducible for a given `seed`, both sequentially and with
  `n_jobs > 1` (worker RNG streams are seeded).
- [`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
  and
  [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
  no longer overwrite the caller’s RNG state.
- Seurat input to the `*_per_sample()` helpers matches images to
  `meta.data` by cell name, so `sample_key` values need not equal the
  image names.
- [`interaction_spot_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/interaction_spot_per_sample.md)
  returns `NA` (not 0 spots) when the spot search fails, and no longer
  requires a `cell` column in `meta.data`.
- Absent cell types give `NA` scores instead of silently dropped rows.

### New features

- **Experimental segmentation-free analysis** of transcript coordinates
  (e.g. Xenium `transcripts.parquet`):
  [`bin_transcripts()`](https://juninamo.github.io/spatialCooccur/reference/bin_transcripts.md),
  model-free cross pair correlation of gene sets
  [`pcf_cross()`](https://juninamo.github.io/spatialCooccur/reference/pcf_cross.md)
  (the relative version is the label-permutation O/E, the continuous
  analogue of `log2_oe`), a random-feature log-Gaussian Cox process
  factor model
  [`fit_spatial_rff()`](https://juninamo.github.io/spatialCooccur/reference/fit_spatial_rff.md)
  after Gundersen, Zhang & Engelhardt (AISTATS 2021) with model-based
  [`rff_pair_correlation()`](https://juninamo.github.io/spatialCooccur/reference/rff_pair_correlation.md),
  and
  [`colocalization_per_sample()`](https://juninamo.github.io/spatialCooccur/reference/colocalization_per_sample.md)
  whose output goes into
  [`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
  with `pair_keys = c("cluster_i", "cluster_j", "r")`. Simulators
  [`simulate_transcripts()`](https://juninamo.github.io/spatialCooccur/reference/simulate_transcripts.md)
  /
  [`simulate_transcripts_groups()`](https://juninamo.github.io/spatialCooccur/reference/simulate_transcripts_groups.md)
  and
  [`lgcp_true_pair_correlation()`](https://juninamo.github.io/spatialCooccur/reference/lgcp_true_pair_correlation.md)
  support validation.
  [`read_xenium_transcripts()`](https://juninamo.github.io/spatialCooccur/reference/read_xenium_transcripts.md)
  reads Xenium transcript tables (binary gene names in older outputs, qv
  filter, gene filtering inside Arrow for 5K panels) and
  [`pcf_matrix()`](https://juninamo.github.io/spatialCooccur/reference/pcf_matrix.md)
  computes all gene-set pairs with cached FFTs.

- [`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md)
  supports paired / repeated-measures designs (e.g. pre- vs
  post-treatment): new `method = "signrank"`, and `method = "perm"` now
  permutes labels within patients when patients appear in both groups
  (the previous between-patient permutation was invalid for such
  designs).

- [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
  also returns `expected` (permutation mean) and `log2_oe` (log2
  observed / expected), an effect size that does not grow with the
  number of cells. Recommended for between-group comparison.

- `*_per_sample()` outputs gain `n_cells`, `n_i`, `n_j` (abundance of
  the two cell types).

- [`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md):
  exact Wilcoxon p-values for small samples; LMM p-values with
  Satterthwaite df (via lmerTest) instead of Wald z; `covariates`;
  `symmetric`; `min_n_per_group`; exact enumeration for small
  permutation tests; a warning on pseudoreplication.

- [`summarize_by_patient()`](https://juninamo.github.io/spatialCooccur/reference/summarize_by_patient.md)
  aggregates image-level results to patients.

- [`generate_sim_groups()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim_groups.md)
  supports several images per patient (`n_images_per_patient`,
  `within_patient_noise`) and group-specific `n_cells`.

- New tutorial: case-control comparison
  (`vignettes/case_control_tutorial.ipynb`).

- Faster counting with sparse matrix products.
