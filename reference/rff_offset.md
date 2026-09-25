# Offset for a residual random-feature model: expression explained by known structure

\*\*Experimental.\*\* Fits, for every gene, a Poisson regression of the
bin counts on covariates that describe structure already known -
cell-type composition of the bins, tissue domains (one-hot), or an
embedding such as PCA, Harmony or SCIGMA - with the bin area as offset,
and returns the fitted log expectation per bin and gene (without the bin
area). Passed as \`offset\` to \[fit_spatial_rff()\], the spatial
factors then capture only the spatially coherent variation that the
covariates do not explain.

## Usage

``` r
rff_offset(binned, covariates, genes = NULL, ridge = 1e-04)
```

## Arguments

- binned:

  A \`binned_transcripts\` object.

- covariates:

  Numeric matrix or data frame with one row per bin of \`binned\` (all
  bins, or only the in-tissue bins). Factors (e.g. domain or section
  labels) are expanded to indicator columns. An intercept is added.

- genes:

  Genes to model (default: all genes of \`binned\`).

- ridge:

  Small ridge penalty stabilising the per-gene fits (collinear or sparse
  covariates).

## Value

A numeric matrix (in-tissue bins x genes, natural log scale) for
\`fit_spatial_rff(offset = )\`.

## See also

\[fit_spatial_rff()\], \[rff_factor_test()\]

## Examples

``` r
tx <- simulate_transcripts(size = 120, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
# known structure: here a smooth coordinate trend as a stand-in covariate
cv <- b$coords[b$coords$in_tissue, c("x", "y")]
off <- rff_offset(b, cv)
fit <- fit_spatial_rff(b, n_factors = 3, offset = off, ard = 2,
                       n_features = 32, max_iter = 50)
fit$factor_strength
#>  factor1  factor2  factor3 
#> 1.357358 1.733956 1.128550 
```
