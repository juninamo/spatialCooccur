# Simulate Spatial Coordinates and Cell Types

Generate synthetic spatial transcriptomics data for simulation and
benchmarking.

## Usage

``` r
generate_sim(
  close_ratio = 0.7,
  n_types = 10,
  max_loc = 800,
  n_perm = 100,
  n_cells = 1500,
  test_type = "circle",
  distance_param = 50,
  seed = 1234
)
```

## Arguments

- close_ratio:

  Proportion of close interactions between selected cell types.

- n_types:

  Number of distinct cell types.

- max_loc:

  Maximum coordinate value (spatial extent).

- n_perm:

  Number of permutations to simulate, for use in future analysis.

- n_cells:

  Total number of cells to simulate.

- test_type:

  Type of spatial pattern to simulate. One of "circle", "line", or
  "distribute".

- distance_param:

  Distance parameter controlling interaction distance.

- seed:

  Random seed for reproducibility.

## Value

A data.frame with simulated spatial coordinates and cell type labels.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
head(df)
#>           x         y   cell_type
#> 1  79.65260 202.11367 cell_type_3
#> 2 111.63717  28.45736 cell_type_1
#> 3 171.85601 147.77884 cell_type_2
#> 4 272.46234 138.46555 cell_type_1
#> 5  60.50458 112.56496 cell_type_1
#> 6 269.51691 297.32977 cell_type_4
table(df$cell_type)
#> 
#> cell_type_1 cell_type_2 cell_type_3 cell_type_4 
#>          76          77          70          77 
```
