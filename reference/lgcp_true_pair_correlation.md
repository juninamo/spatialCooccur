# True cross pair correlation of a simulated LGCP

Closed-form cross pair correlation between two gene sets for data from
\[simulate_transcripts()\], for validating \[pcf_cross()\] and
\[rff_pair_correlation()\]. See \[rff_pair_correlation()\] for the three
types.

## Usage

``` r
lgcp_true_pair_correlation(
  truth,
  set_a,
  set_b,
  r,
  type = c("relative", "full", "composition")
)
```

## Arguments

- truth:

  The \`truth\` attribute of a \[simulate_transcripts()\] result.

- set_a, set_b:

  Gene-set names.

- r:

  Distances (micrometres).

- type:

  \`"relative"\`, \`"full"\` or \`"composition"\`.

## Value

A data.frame with \`r\`, \`g\` and \`log_g\`.

## Examples

``` r
tx <- simulate_transcripts(size = 100, rate = 0.01)
lgcp_true_pair_correlation(attr(tx, "truth"), "A", "B", r = c(5, 10, 20))
#>    r         g       log_g
#> 1  5 0.7870381 -0.23947867
#> 2 10 0.8471166 -0.16591695
#> 3 20 0.9874068 -0.01267321
```
