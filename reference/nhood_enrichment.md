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

A list with matrices \`zscore\`, \`count\` (observed), \`expected\`
(mean of the permutation null) and \`log2_oe\` (log2 observed /
expected, with a pseudocount of one mean edge weight). Unlike the
z-score, which grows with the number of cells, \`log2_oe\` is an effect
size that is comparable across samples of different size.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
res <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = 10,
                        n_perms = 50, n_jobs = 1)
round(res$zscore, 1)
#>                    Clustercell_type_1 Clustercell_type_2 Clustercell_type_3
#> Clustercell_type_1               -0.9                1.6                1.2
#> Clustercell_type_2                1.1                0.4               -2.5
#> Clustercell_type_3                1.2               -1.9                1.4
#> Clustercell_type_4               -0.7                0.3                0.1
#>                    Clustercell_type_4
#> Clustercell_type_1               -1.9
#> Clustercell_type_2               -1.3
#> Clustercell_type_3               -0.8
#> Clustercell_type_4                2.5
round(res$log2_oe, 2)
#>                    Clustercell_type_1 Clustercell_type_2 Clustercell_type_3
#> Clustercell_type_1              -0.08               0.15               0.11
#> Clustercell_type_2               0.10               0.03              -0.26
#> Clustercell_type_3               0.12              -0.21               0.12
#> Clustercell_type_4              -0.09               0.03               0.01
#>                    Clustercell_type_4
#> Clustercell_type_1              -0.26
#> Clustercell_type_2              -0.14
#> Clustercell_type_3              -0.09
#> Clustercell_type_4               0.21
```
