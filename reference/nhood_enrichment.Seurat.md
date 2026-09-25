# Neighborhood Enrichment (Seurat Method)

Neighborhood Enrichment (Seurat Method)

## Usage

``` r
nhood_enrichment.Seurat(
  seurat_obj,
  cluster_key,
  neighbors.k = 30,
  connectivity_key = "nn",
  transformation = FALSE,
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

  If \`TRUE\`, each link from cell u is weighted 1 / (1 + d_u), with d_u
  the number of cells that chose u as a neighbour (a mild down-weighting
  of hub cells). \`FALSE\` (default since 0.99.3) counts every kNN link
  once, as squidpy does. With a kNN graph every cell sends k links at
  any density, so the weighting is not a density correction (density is
  handled by the label shuffles); in simulations it kept the calibration
  but lowered power.

- n_perms:

  Number of permutations for significance testing.

- seed:

  Random seed for reproducibility.

- n_jobs:

  Number of cores to use in parallel.

## Value

Updated Seurat object; \`misc\[\[paste0(cluster_key,
"\_nhood_enrichment")\]\]\` holds \`zscore\`, \`count\`, \`expected\`
(permutation mean), \`log2_oe\`, \`log2_oe_raw\`, \`pvalue\`, \`padj\`
and \`padj_bh\` (see \[nhood_enrichment()\]).

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
#> Clustercell_type_1               -0.5                2.2                1.8
#> Clustercell_type_2                1.4                2.7               -2.3
#> Clustercell_type_3                1.8               -1.2                1.1
#> Clustercell_type_4               -1.7               -0.9                0.1
#>                    Clustercell_type_4
#> Clustercell_type_1               -3.5
#> Clustercell_type_2               -3.3
#> Clustercell_type_3               -1.2
#> Clustercell_type_4                2.4
```
