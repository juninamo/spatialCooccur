# Simulate a multi-sample case-control transcript dataset

\*\*Experimental.\*\* Calls \[simulate_transcripts()\] for every sample,
with a group-specific niche loading for the co-localizing gene sets plus
patient-level noise.

## Usage

``` r
simulate_transcripts_groups(
  n_samples_per_group = 4,
  group_coloc = list(control = 0.3, case = 0.9),
  coloc_sets = c("A", "B"),
  between_sample_sd = 0.1,
  seed = 1,
  ...
)
```

## Arguments

- n_samples_per_group:

  Number of samples (patients) per group.

- group_coloc:

  Named list, one numeric niche loading per group, applied to every set
  in \`coloc_sets\`, e.g. \`list(control = 0.3, case = 0.9)\`.

- coloc_sets:

  Gene sets that share the niche field.

- between_sample_sd:

  SD of the per-sample noise on the niche loading.

- seed:

  Random seed.

- ...:

  Passed to \[simulate_transcripts()\].

## Value

A data.frame as \[simulate_transcripts()\] with \`sample_id\`, \`group\`
and \`patient\` columns.

## Examples

``` r
tx <- simulate_transcripts_groups(n_samples_per_group = 2, size = 120,
                                  rate = 0.01)
table(tx$sample_id, tx$group)
#>            
#>             case control
#>   case_1    3387       0
#>   case_2    2165       0
#>   control_1    0    5409
#>   control_2    0    9245
```
