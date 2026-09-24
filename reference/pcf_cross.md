# Empirical cross pair correlation of two gene sets (segmentation-free)

\*\*Experimental.\*\* Estimate the cross pair correlation function
\\g\_{AB}(r)\\ between the transcripts of gene set A and gene set B from
binned counts, with FFT-based edge correction over the in-tissue bins.
\\g\_{AB}(r) \> 1\\ means B transcripts are found at distance \`r\` from
A transcripts more often than under independence (the continuous
analogue of the neighborhood enrichment \`log2_oe\`).

## Usage

``` r
pcf_cross(
  binned,
  set_a,
  set_b = set_a,
  r_max = 100,
  r_step = binned$grid$bin_size,
  relative = FALSE
)
```

## Arguments

- binned:

  A \`binned_transcripts\` object from \[bin_transcripts()\].

- set_a, set_b:

  Character vectors of genes (B defaults to A, giving the auto pair
  correlation).

- r_max:

  Largest distance to report.

- r_step:

  Width of the distance annuli (default: the bin size).

- relative:

  If \`TRUE\`, divide by the pair correlation of all transcripts. This
  is the observed / expected ratio of A-B pairs at distance \`r\` when
  gene labels are permuted over the fixed transcript positions (the
  normalized mark connection function) - the segmentation-free
  counterpart of the label-permutation \`log2_oe\` of
  \[nhood_enrichment()\]. Variation in overall cellularity cancels in
  the ratio, so it is the recommended value for comparing samples.

## Value

A data.frame with \`r\` (annulus midpoint), \`g\`, \`log_g\` and
\`n_lags\` (number of grid lags averaged).

## Examples

``` r
tx <- simulate_transcripts(size = 150, rate = 0.01, coloc = c(A = 1, B = 1))
b <- bin_transcripts(tx, bin_size = 3)
genes_a <- unique(tx$gene[tx$gene_set == "A"])
genes_b <- unique(tx$gene[tx$gene_set == "B"])
head(pcf_cross(b, genes_a, genes_b, r_max = 30))
#>      r        g      log_g n_lags
#> 1  1.5 1.677433 0.51726487      1
#> 2  4.5 1.550645 0.43867105      8
#> 3  7.5 1.408187 0.34230288     16
#> 4 10.5 1.249410 0.22267179     20
#> 5 13.5 1.159685 0.14814809     24
#> 6 16.5 1.040379 0.03958543     40
```
