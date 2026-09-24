# Per-sample segmentation-free co-localization of gene sets

\*\*Experimental.\*\* For every sample, bin the transcripts and compute
the log cross pair correlation between gene sets at the requested
distances, either model-free (\[pcf_cross()\]) or from a fitted
random-feature LGCP (\[fit_spatial_rff()\] +
\[rff_pair_correlation()\]). The tidy output goes straight into
\[compare_groups()\] with \`value = "log_g"\` and \`pair_keys =
c("cluster_i", "cluster_j", "r")\`.

## Usage

``` r
colocalization_per_sample(
  tx,
  sample_key,
  group_key,
  patient_key = NULL,
  gene_sets,
  pairs = NULL,
  r = c(10, 25, 50),
  method = c("empirical", "model"),
  bin_size = 4,
  relative = TRUE,
  model_type = c("composition", "relative", "full"),
  n_cores = 1,
  ...
)
```

## Arguments

- tx:

  Transcript table with a sample column.

- sample_key, group_key, patient_key:

  Column names.

- gene_sets:

  Named list of gene vectors (e.g. cell-type markers).

- pairs:

  Two-column data.frame or list of length-2 character vectors of
  gene-set names to test (default: all unordered pairs, including
  self-pairs).

- r:

  Distances at which to report \`log_g\`.

- method:

  \`"empirical"\` or \`"model"\`.

- bin_size:

  Bin size passed to \[bin_transcripts()\].

- relative:

  For \`method = "empirical"\`: if \`TRUE\` (default), report the
  relative pair correlation (label-permutation O/E; see
  \[pcf_cross()\]), which is robust to differences in overall
  cellularity between samples.

- model_type:

  For \`method = "model"\`: the \`type\` passed to
  \[rff_pair_correlation()\] (default \`"composition"\`).

- n_cores:

  Number of samples processed in parallel (forked processes via
  \[parallel::mclapply()\]; sequential on Windows).

- ...:

  Passed to \[fit_spatial_rff()\] (for \`method = "model"\`) or
  \[bin_transcripts()\].

## Value

A data.frame with \`sample_id\`, \`cluster_i\`, \`cluster_j\`, \`r\`,
\`log_g\`, \`n_transcripts\`, \`group\`, \`patient\`.

## Examples

``` r
tx <- simulate_transcripts_groups(n_samples_per_group = 2, size = 120,
                                  rate = 0.01)
sets <- split(unique(tx[, c("gene", "gene_set")])$gene,
              unique(tx[, c("gene", "gene_set")])$gene_set)
res <- colocalization_per_sample(tx, "sample_id", "group", "patient",
                                 gene_sets = sets,
                                 pairs = list(c("A", "B")), r = c(10, 20))
res
#>   sample_id cluster_i cluster_j  r       log_g n_transcripts   group   patient
#> 1 control_1         A         B 10 -0.12979684          5409 control control_1
#> 2 control_1         A         B 20 -0.02837551          5409 control control_1
#> 3 control_2         A         B 10 -0.37804801          9245 control control_2
#> 4 control_2         A         B 20 -0.19748333          9245 control control_2
#> 5    case_1         A         B 10 -0.08906563          3387    case    case_1
#> 6    case_1         A         B 20 -0.01324909          3387    case    case_1
#> 7    case_2         A         B 10 -0.41569010          2165    case    case_2
#> 8    case_2         A         B 20 -0.25463029          2165    case    case_2
```
