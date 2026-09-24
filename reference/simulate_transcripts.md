# Simulate transcript coordinates from a multivariate log-Gaussian Cox process

\*\*Experimental.\*\* Genes are organized in gene sets (e.g. marker
genes of cell types). The log-intensity of gene \`j\` in set \`S\` at
location \`u\` is \$\$\log\lambda_j(u) = \log\rho + \sigma_d d(u) +
\sigma_s f_S(u) + c_S n(u),\$\$ where \`d\` is a shared cellularity
field, \`f_S\` a set-specific territory field and \`n\` a "niche" field
shared by the sets in \`coloc\`. All fields are unit-variance Gaussian
random fields with an RBF covariance, generated exactly by FFT (not with
random features). Transcripts are Poisson counts on a \`resolution\`
grid, placed uniformly within each grid cell.

## Usage

``` r
simulate_transcripts(
  gene_sets = c("A", "B", "C"),
  n_genes_per_set = 5,
  size = 500,
  resolution = 2,
  rate = 0.02,
  density_sd = 0.5,
  density_lengthscale = 80,
  set_sd = 1,
  set_lengthscale = 15,
  coloc = c(A = 0.8, B = 0.8),
  niche_lengthscale = 20,
  seed = 1
)
```

## Arguments

- gene_sets:

  Character vector of gene-set names.

- n_genes_per_set:

  Number of genes in each set.

- size:

  Side length of the square window (micrometres).

- resolution:

  Simulation grid spacing (micrometres).

- rate:

  Baseline transcripts per square micrometre per gene.

- density_sd, density_lengthscale:

  SD and length scale of the shared cellularity field.

- set_sd, set_lengthscale:

  SD and length scale of the set-specific fields.

- coloc:

  Named numeric vector of niche loadings, e.g. \`c(A = 0.8, B = 0.8)\`
  makes sets A and B co-localize; opposite signs make them segregate.

- niche_lengthscale:

  Length scale of the niche field.

- seed:

  Random seed.

## Value

A data.frame with columns \`x\`, \`y\`, \`gene\`, \`gene_set\`, carrying
an attribute \`truth\` (the generating parameters) usable with
\[lgcp_true_pair_correlation()\].

## Examples

``` r
tx <- simulate_transcripts(size = 150, rate = 0.01, coloc = c(A = 1, B = 1))
head(tx)
#>           x         y gene gene_set
#> 1  64.60226 0.1340603  A_1        A
#> 2  87.86456 0.2502116  A_1        A
#> 3  89.92663 0.5943523  A_1        A
#> 4  88.75204 0.9862053  A_1        A
#> 5 101.01818 0.4752379  A_1        A
#> 6 131.72347 0.8982624  A_1        A
table(tx$gene_set)
#> 
#>    A    B    C 
#> 1956 1642 1111 
```
