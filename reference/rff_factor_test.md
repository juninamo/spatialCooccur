# Significance of residual spatial factors by parametric bootstrap

\*\*Experimental.\*\* Tests whether the factors of a
\[fit_spatial_rff()\] fit are stronger than factors fitted to data
without any factor structure. Counts are simulated from the fitted null
model (bin area, offset, gene intercepts, cellularity field and
dispersion, but no factors), the same model is refitted to each
simulated data set, and the largest factor strength (norm of a factor's
loadings) is recorded. Each observed factor is compared with this
max-null distribution, which controls the family-wise error over
factors.

## Usage

``` r
rff_factor_test(fit, binned, n_boot = 19, max_iter = NULL, seed = 1)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\].

- binned:

  The \`binned_transcripts\` object used for the fit.

- n_boot:

  Number of simulated null data sets.

- max_iter:

  Iterations for the refits (default: those of the fit).

- seed:

  Random seed.

## Value

A data frame with one row per factor: \`factor\`, \`strength\`,
\`lengthscale\`, \`p\` (share of null data sets whose strongest factor
is at least as strong, with the +1 correction), and the top genes by
absolute loading. The null maxima are attached as attribute
\`null_max\`.

## See also

\[fit_spatial_rff()\], \[rff_offset()\], \[rff_fields()\]

## Examples

``` r
# \donttest{
tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
fit <- fit_spatial_rff(b, n_factors = 2, ard = 2, n_features = 24, max_iter = 40)
rff_factor_test(fit, b, n_boot = 4)
#>    factor strength lengthscale   p
#> 1 factor1 1.773590    10.17890 0.2
#> 2 factor2 1.479139    21.97312 0.2
#>                                                         top_genes
#> 1 A_3 (-0.71), B_3 (-0.70), A_1 (-0.66), A_2 (-0.65), B_1 (-0.63)
#> 2 B_2 (+0.94), B_3 (+0.82), B_1 (+0.69), C_1 (+0.20), A_1 (-0.19)
# }
```
