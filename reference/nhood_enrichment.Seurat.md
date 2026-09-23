# Neighborhood Enrichment (Seurat Method)

Neighborhood Enrichment (Seurat Method)

## Usage

``` r
nhood_enrichment.Seurat(
  seurat_obj,
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

- seurat_obj:

  A Seurat object with spatial coordinates.

- cluster_key:

  Metadata column for cluster IDs.

- neighbors.k:

  Number of neighbors to construct graph.

- connectivity_key:

  Which graph to use: "nn" or "snn".

- transformation:

  Logical, whether to normalize adjacency matrix.

- n_perms:

  Number of permutations for significance testing.

- seed:

  Random seed for reproducibility.

- n_jobs:

  Number of cores to use in parallel.

## Value

Updated Seurat object; \`misc\[\[paste0(cluster_key,
"\_nhood_enrichment")\]\]\` holds \`zscore\`, \`count\`, \`expected\`
(permutation mean) and \`log2_oe\`.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
df$sample_id <- "fov1"
seu <- sim_to_seurat(df)
seu <- nhood_enrichment.Seurat(seu, cluster_key = "cell_type",
                               neighbors.k = 10, n_perms = 20, n_jobs = 1)
res <- SeuratObject::Misc(seu, slot = "cell_type_nhood_enrichment")
round(res$zscore, 1)
#>                    Clustercell_type_1 Clustercell_type_2 Clustercell_type_3
#> Clustercell_type_1               -1.0                1.7                1.1
#> Clustercell_type_2                1.1                0.2               -2.5
#> Clustercell_type_3                1.2               -2.1                1.0
#> Clustercell_type_4               -1.0                0.5                0.2
#>                    Clustercell_type_4
#> Clustercell_type_1               -2.7
#> Clustercell_type_2               -1.4
#> Clustercell_type_3               -0.8
#> Clustercell_type_4                2.6
```
