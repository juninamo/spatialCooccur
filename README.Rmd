
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialCooccur <img src="man/figures/logo.png" align="right" height="138" />

`spatialCooccur` is an R package for analyzing spatial co-occurrence and
neighborhood interactions in spatial transcriptomics data. It is built
around Seurat objects and provides tools to compute co-occurrence
enrichment, perform permutation-based tests, and visualize local
interaction scores.

## Installation

You can install the development version from GitHub using:

``` r
# install.packages("devtools")
devtools::install_github("juninamo/spatialCooccur")
```

see details
[here](https://github.com/juninamo/Spatial_Neighborhood_Analysis)

## Features

- Simulate spatial transcriptomic layouts with `generate_sim()`
- Calculate neighborhood co-occurrence enrichment with
  `nhood_enrichment()`
- Identify local interaction zones using `cooccur_local()`
- Permutation-based z-score estimation
- Works seamlessly with Seurat spatial objects

## Example

``` r
library(spatialCooccur)
library(Seurat)

# Simulate example data
df <- generate_sim(n_cells = 500, n_types = 4, test_type = "circle")

# Compute neighborhood enrichment
res <- nhood_enrichment(df, cluster_key = "cell_type")

# View z-score matrix
res$zscore
```

## License

MIT © Jun Inamo
