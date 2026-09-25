# Compute Co-occurrence Count Matrix

Compute Co-occurrence Count Matrix

## Usage

``` r
compute_count(
  adj,
  int_clust_row,
  int_clust_col,
  n_cls,
  cluster_data,
  transformation = FALSE
)
```

## Arguments

- adj:

  Adjacency matrix.

- int_clust_row:

  Vector of cluster labels for rows.

- int_clust_col:

  Vector of cluster labels for columns.

- n_cls:

  Number of clusters.

- cluster_data:

  Original cluster assignments.

- transformation:

  If \`TRUE\`, \`adj\` holds link weights (e.g. 1 / (1 + d_u)) that are
  summed; \`FALSE\` (default) counts binary links.

## Value

A co-occurrence count matrix.

## Examples

``` r
set.seed(1)
adj <- Matrix::rsparsematrix(20, 20, density = 0.2)
cl <- factor(sample(c("a", "b"), 20, replace = TRUE))
lab <- paste0("Cluster", cl)
compute_count(adj, lab, lab, n_cls = 2, cluster_data = cl)
#>          Clustera Clusterb
#> Clustera        0        1
#> Clusterb        0        1
```
