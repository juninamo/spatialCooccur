# Compute Enrichment Ratios from Count Matrix

Compute Enrichment Ratios from Count Matrix

## Usage

``` r
compute_co_occurrence_ratio(co_occur_count)
```

## Arguments

- co_occur_count:

  Matrix of observed co-occurrence counts.

## Value

A matrix of normalized enrichment ratios.

## Examples

``` r
counts <- matrix(c(10, 2, 2, 6), 2,
                 dimnames = list(c("A", "B"), c("A", "B")))
compute_co_occurrence_ratio(counts)
#>           A         B
#> A 1.3888889 0.4166667
#> B 0.4166667 1.8750000
```
