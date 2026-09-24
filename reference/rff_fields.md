# Latent fields of a fitted random-feature LGCP at cells, transcripts or any point

\*\*Experimental.\*\* \[fit_spatial_rff()\] represents the cellularity
field and the \`K\` spatial factors as smooth functions of the
coordinates (random Fourier features), so they can be evaluated
anywhere, not only on the bin grid stored in \`fit\$field_grid\`.
\`rff_fields()\` evaluates them at the given coordinates - cell
centroids, transcripts or any other points - and optionally averages the
values per group (e.g. per cell, from the transcripts assigned to it).

## Usage

``` r
rff_fields(
  fit,
  coords,
  x_col = "x",
  y_col = "y",
  by = NULL,
  which = c("all", "factors", "density"),
  warn_outside = TRUE
)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\].

- coords:

  A data frame (or matrix) with the coordinates, in the same unit and
  frame as the transcripts used for the fit: one row per cell centroid,
  transcript or point.

- x_col, y_col:

  Column names of the coordinates in \`coords\`.

- by:

  Optional column of \`coords\` (e.g. \`"cell_id"\`) whose values define
  groups; the fields are then averaged within each group, e.g. over the
  transcripts of a cell. \`NULL\` (default) returns one row per point.

- which:

  \`"all"\` (default), \`"factors"\` or \`"density"\`.

- warn_outside:

  Warn when points lie outside the fitted tissue (farther than one bin
  from any in-tissue bin), where the fields are extrapolated.

## Value

A data frame with one row per point (or per group when \`by\` is given,
with the group in the first column and the number of points in
\`n_points\`), and one column per field: \`density\` (if fitted) and
\`factor1\` ... \`factorK\`. Row names are kept from \`coords\` when
\`by\` is \`NULL\`.

## Details

The fields are on the unit-variance scale of \`fit\$field_grid\` (each
field has variance 1 over the fitted tissue). The contribution of the
factors to the log intensity of gene \`j\` at a point is \`fields\[,
factors\] (plus \`fit\$sigma0\` times the density field when it was
fitted).

## See also

\[fit_spatial_rff()\], \[rff_pair_correlation()\]

## Examples

``` r
tx <- simulate_transcripts(size = 120, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
fit <- fit_spatial_rff(b, n_factors = 2, n_features = 32, max_iter = 50)
# one row per transcript
head(rff_fields(fit, tx))
#>     density    factor1    factor2
#> 1 0.4901363 -0.3997969 -0.9228010
#> 2 0.5209908 -0.3979304 -0.8197534
#> 3 0.5966679 -0.4585247 -0.6593935
#> 4 0.5924033 -0.4549208 -0.6671474
#> 5 0.6576537 -0.3733627 -0.3402887
#> 6 0.7173196 -0.3899766 -0.2344311
# averaged per "cell" (here: 20 x 20 um squares as stand-in cells)
tx$cell_id <- paste(floor(tx$x / 20), floor(tx$y / 20))
head(rff_fields(fit, tx, by = "cell_id"))
#>   cell_id n_points   density     factor1     factor2
#> 1     0 0      335 0.8294371 -0.49487114 -0.06366131
#> 2     0 1      429 0.8669994 -0.87635103  0.38787658
#> 3     0 2      194 0.6624920 -0.82060599  0.58997225
#> 4     0 3       73 0.4368530 -0.08299567  1.44415322
#> 5     0 4       92 0.3883809 -0.51402333  1.39180678
#> 6     0 5      126 0.3127352 -0.38058497  0.70779222
```
