# Associate co-localization with a continuous clinical variable

For every cell-type pair (or pair x distance), tests whether a
co-localization score from a \`\*\_per_sample()\` helper changes with a
continuous variable such as CRP, a disease-activity score or age. Images
of the same patient are averaged (\`"spearman"\`, \`"lm"\`, \`"perm"\`)
or modelled with a random patient intercept (\`"lmm"\`), so the patient
is the unit of analysis.

## Usage

``` r
associate_continuous(
  per_sample_df,
  variable,
  value = "log2_oe",
  method = c("spearman", "lm", "lmm", "perm"),
  patient_key = "patient",
  covariates = NULL,
  pair_keys = NULL,
  n_perms = 2000,
  adjust = "BH",
  min_n = 4,
  seed = 1234
)
```

## Arguments

- per_sample_df:

  Output of a \`\*\_per_sample()\` helper (one row per sample x pair),
  with a \`patient\` column when patients have several images.

- variable:

  Name of the numeric column with the clinical variable, or a named
  numeric vector indexed by patient (or sample) ID.

- value:

  Score column to test, e.g. \`"log2_oe"\` or \`"log_g_rel"\`.

- method:

  One of \`"spearman"\`, \`"lm"\`, \`"lmm"\`, \`"perm"\`.

- patient_key:

  Column identifying patients (default \`"patient"\`).

- covariates:

  Optional covariate columns for \`"lm"\` / \`"lmm"\`.

- pair_keys:

  Columns identifying one test (default \`cluster_i\`, \`cluster_j\`,
  plus \`r\` when present).

- n_perms:

  Permutations for \`method = "perm"\`.

- adjust:

  Multiple-testing correction passed to \[stats::p.adjust()\].

- min_n:

  Minimum number of patients with both values.

- seed:

  Random seed for \`"perm"\`.

## Value

A data.frame with one row per pair: \`n_patients\`, \`estimate\`
(Spearman rho, or the slope per unit of \`variable\`), \`p\`, \`padj\`,
\`method\`.

## Details

\* \`"spearman"\`: Spearman correlation between patient means and the
variable (exact p-value for small cohorts without ties). Robust to
outliers and to non-linear but monotone relationships. \* \`"lm"\`:
linear regression of patient means on the variable, optionally adjusted
for \`covariates\`; reports the slope. \* \`"lmm"\`: linear mixed model
on images, \`score ~ variable + covariates + (1 \| patient)\` (lmerTest,
Satterthwaite df); uses every image when patients have several. \*
\`"perm"\`: Spearman correlation of patient means with a permutation
p-value (the variable is shuffled between patients).

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 8,
                          group_close_ratio = list(all = 0.4),
                          n_types = 3, n_cells = 300, max_loc = 300,
                          test_type = "distribute", distance_param = 10,
                          seed = 2)
ps <- nhood_enrichment_per_sample(df, sample_key = "sample_id",
                                  group_key = "group",
                                  cluster_key = "cell_type",
                                  patient_key = "patient",
                                  neighbors.k = 10, n_perms = 30, n_jobs = 1)
crp <- setNames(rexp(8), unique(ps$patient))
head(associate_continuous(ps, crp, value = "log2_oe"))
#>     cluster_i   cluster_j n_patients   estimate         p      padj   method
#> 6 cell_type_3 cell_type_2          8  0.5000000 0.2161706 0.7959821 spearman
#> 8 cell_type_2 cell_type_3          8  0.4761905 0.2430556 0.7959821 spearman
#> 2 cell_type_2 cell_type_1          8 -0.4285714 0.2992063 0.7959821 spearman
#> 4 cell_type_1 cell_type_2          8 -0.3333333 0.4278770 0.7959821 spearman
#> 5 cell_type_2 cell_type_2          8 -0.2380952 0.5821429 0.7959821 spearman
#> 9 cell_type_3 cell_type_3          8 -0.2380952 0.5821429 0.7959821 spearman
```
