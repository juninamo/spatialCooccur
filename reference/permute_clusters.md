# Permute Cluster Assignments and Recompute Counts

Permute Cluster Assignments and Recompute Counts

## Usage

``` r
permute_clusters(adj, int_clust, n_cls, cluster_data, transformation)
```

## Arguments

- adj:

  Adjacency matrix.

- int_clust:

  Cluster labels.

- n_cls:

  Number of clusters.

- cluster_data:

  Original cluster data.

- transformation:

  Whether to apply adjacency transformation.

  Cluster labels are shuffled once and the same permutation is applied
  to rows and columns of the adjacency matrix.

## Value

Permuted co-occurrence count matrix.

## Examples

``` r
set.seed(1)
adj <- Matrix::rsparsematrix(20, 20, density = 0.2)
cl <- factor(sample(c("a", "b"), 20, replace = TRUE))
permute_clusters(adj, paste0("Cluster", cl), n_cls = 2, cluster_data = cl,
                 transformation = TRUE)
#>          Clustera Clusterb
#> Clustera    0.260   5.1600
#> Clusterb   -1.889  -6.7143
```
