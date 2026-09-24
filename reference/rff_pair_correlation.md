# Model-based cross pair correlation from a fitted random-feature LGCP

\*\*Experimental.\*\* Under an LGCP the pair correlation of genes \`j\`
and \`j'\` is \\g\_{jj'}(h) = \exp(C\_{jj'}(h))\\, with \\C\_{jj'}\\ the
cross-covariance of their log-intensities. From the fit, \\C\_{jj'}(h) =
\tilde L_j \hat\Sigma(h) \tilde L\_{j'}^\top\\, where \\\tilde L\\ are
the loadings (with \\\sigma_0\\ for the cellularity field) and
\\\hat\Sigma(h)\\ is the empirical cross-covariance matrix of the fitted
unit-variance fields at distance \`h\` (plug-in estimate). Using the
fields' actual cross-covariances, rather than assuming independent
factors, keeps the estimate correct when fitted factors are correlated.
For gene sets, gene pairs are weighted by their expected intensities.

## Usage

``` r
rff_pair_correlation(
  fit,
  set_a,
  set_b,
  r,
  type = c("relative", "full", "composition"),
  method = c("intensity", "gaussian")
)
```

## Arguments

- fit:

  A \`spatial_rff_fit\` from \[fit_spatial_rff()\].

- set_a, set_b:

  Character vectors of genes.

- r:

  Distances.

- type:

  \`"relative"\` (default), \`"full"\` or \`"composition"\`.

- method:

  \`"intensity"\` (default): apply the FFT pair-correlation estimator of
  \[pcf_cross()\] to the fitted intensity surfaces
  \\\hat\lambda_j(u)\\ - a denoised version of the empirical estimate
  that makes no Gaussian assumption. \`"gaussian"\`: the closed-form
  LGCP value \\\exp(C\_{jj'}(h))\\, which is exact for Gaussian
  log-intensity fields but can be far off for sparse, strongly clustered
  genes (vessels, meninges) in real tissue.

## Value

A data.frame with \`r\`, \`g\`, \`log_g\`.

## Details

Three summaries are available: \* \`"relative"\`: \\g\_{AB}(h) /
g\_{TT}(h)\\, with T all modelled transcripts - the model-based
counterpart of \`pcf_cross(relative = TRUE)\` and of the
label-permutation \`log2_oe\`. \* \`"full"\`: \\g\_{AB}(h)\\, including
the shared cellularity field. \* \`"composition"\`: \\g\_{AB}(h)\\
without the cellularity field.

\*\*Which summary to compare between groups.\*\* In simulations (6 vs 6
samples, Wilcoxon test): \* \`pcf_cross(relative = TRUE)\` (model-free)
detected a planted change and stayed calibrated when groups differed
only in cellularity heterogeneity; \`"full"\` gave false positives. \*
\`type = "composition"\` with the default \`method = "intensity"\` and
\`offset = "area"\` was calibrated and the most powerful when the panel
contained many cell types (8 gene sets). With only 3 gene sets, two of
which co-localize, the shared niche is hard to tell from the cellularity
field and occasional false positives remained; there, fit with \`offset
= "smoothed_total"\` and \`offset_genes\` set to reference genes outside
the tested pair (calibrated, but less powerful). \* \`method =
"intensity"\` estimates are shrunk towards 0 at short distances
(conservative); \`method = "gaussian"\` is unbiased when the
log-intensities are Gaussian but was badly miscalibrated on real Xenium
data.

## Examples

``` r
tx <- simulate_transcripts(size = 120, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
fit <- fit_spatial_rff(b, n_factors = 2, n_features = 32, max_iter = 50)
rff_pair_correlation(fit, c("A_1", "A_2"), c("B_1", "B_2"), r = c(5, 10, 20))
#>    r         g        log_g
#> 1  5 1.0044906  0.004480519
#> 2 10 0.9909704 -0.009070575
#> 3 20 0.9336334 -0.068671429
```
