# Cross pair correlation for every pair of gene sets

\*\*Experimental.\*\* Compute \[pcf_cross()\] for all pairs of gene sets
at once. The FFT of each set is computed a single time, which makes
whole 5K-panel Xenium sections with dozens of pairs practical.

## Usage

``` r
pcf_matrix(
  binned,
  gene_sets,
  r_max = 100,
  r_step = binned$grid$bin_size,
  reference = binned$genes
)
```

## Arguments

- binned:

  A \`binned_transcripts\` object.

- gene_sets:

  Named list of gene vectors.

- r_max, r_step:

  Distance range and annulus width.

- reference:

  Genes defining the reference pattern for the relative pair correlation
  (default: all binned genes; with a marker-only binning this is "all
  marker transcripts", the analogue of permuting labels among typed
  cells).

## Value

A data.frame with \`cluster_i\`, \`cluster_j\` (unordered, i \<= j in
the order of \`gene_sets\`), \`r\`, \`log_g\` (full) and \`log_g_rel\`
(relative to the reference).

## Examples

``` r
tx <- simulate_transcripts(size = 150, rate = 0.01)
b <- bin_transcripts(tx, bin_size = 5)
sets <- split(attr(tx, "truth")$genes, attr(tx, "truth")$set_of)
head(pcf_matrix(b, sets, r_max = 30))
#>   cluster_i cluster_j    r       log_g   log_g_rel  n_i  n_j
#> 1         A         A  2.5  0.85805868  0.36733670 1809 1809
#> 2         A         A  7.5  0.72673726  0.30208432 1809 1809
#> 3         A         A 12.5  0.49567244  0.19493390 1809 1809
#> 4         A         A 17.5  0.23823740  0.06251405 1809 1809
#> 5         A         A 22.5 -0.02694136 -0.07734969 1809 1809
#> 6         A         A 27.5 -0.28189769 -0.19444367 1809 1809
```
