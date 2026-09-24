# Neighborhood Enrichment (Generic method)

Neighborhood Enrichment (Generic method)

## Usage

``` r
nhood_enrichment(
  df,
  cluster_key,
  neighbors.k = 30,
  connectivity_key = "nn",
  transformation = TRUE,
  n_perms = 100,
  seed = 1938493,
  n_jobs = 4
)
```

## Arguments

- df:

  Data.frame with spatial and cluster metadata.

- cluster_key:

  Column with cluster labels.

- neighbors.k:

  Number of neighbors to use.

- connectivity_key:

  Type of graph: "nn" or "snn".

- transformation:

  Whether to normalize adjacency matrix.

- n_perms:

  Number of permutations.

- seed:

  Random seed.

- n_jobs:

  Number of parallel jobs. \`1\` runs sequentially.

## Value

A list of cell-type x cell-type matrices: \* \`log2_oe\`: effect size,
log2 observed / expected, centred on the label shuffles so that it is 0
on average without interaction for any number of cells (the log of a
ratio of small counts is otherwise biased slightly below 0 for rare cell
types). Unlike the z-score, which grows with the number of cells, it is
comparable across samples. \* \`log2_oe_raw\`: log2((count + c) /
(expected + c)) without centring, with c one mean edge weight. \*
\`pvalue\`: within-sample test per unordered pair (contacts i -\> j and
j -\> i summed), two-sided normal p-value from the shuffles; use it for
a single pre-specified pair. \* \`padj\`: family-wise adjusted p-value
over all K (K + 1) / 2 pairs by the Westfall-Young max-T permutation
method (the share of shuffles whose largest \|z\| reaches the observed
one). Calibrated for any number of cell types, including rare ones; its
smallest value is 1 / (n_perms + 1). \* \`padj_bh\`: Benjamini-Hochberg
on \`pvalue\` (can exceed its level when there are many rare cell
types). \* \`zscore\`, \`count\` (observed), \`expected\` (mean of the
shuffles).

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
res <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = 10,
                        n_perms = 50, n_jobs = 1)
round(res$zscore, 1)
#>                    Clustercell_type_1 Clustercell_type_2 Clustercell_type_3
#> Clustercell_type_1               -0.5                1.0                1.8
#> Clustercell_type_2                0.4                1.3               -2.8
#> Clustercell_type_3                1.0               -2.2                1.5
#> Clustercell_type_4               -0.5               -0.1                0.5
#>                    Clustercell_type_4
#> Clustercell_type_1               -2.0
#> Clustercell_type_2               -2.4
#> Clustercell_type_3               -0.9
#> Clustercell_type_4                3.0
round(res$log2_oe, 2)
#>                    Clustercell_type_1 Clustercell_type_2 Clustercell_type_3
#> Clustercell_type_1              -0.04               0.11               0.16
#> Clustercell_type_2               0.04               0.12              -0.30
#> Clustercell_type_3               0.10              -0.25               0.13
#> Clustercell_type_4              -0.05              -0.01               0.07
#>                    Clustercell_type_4
#> Clustercell_type_1              -0.25
#> Clustercell_type_2              -0.27
#> Clustercell_type_3              -0.09
#> Clustercell_type_4               0.25
```
