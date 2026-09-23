# Compare disease groups across samples

Given a tidy per-sample data.frame produced by one of the
\`\*\_per_sample()\` helpers, run a per-cluster-pair statistical test
between groups and return tidy results with multiple-testing adjustment.

## Usage

``` r
compare_groups(
  per_sample_df,
  value = "zscore",
  group_key = "group",
  patient_key = NULL,
  method = c("wilcox", "t", "lmm", "perm"),
  n_perms = 1000,
  adjust = "BH",
  pair_keys = c("cluster_i", "cluster_j"),
  ref_group = NULL,
  covariates = NULL,
  symmetric = FALSE,
  min_n_per_group = 2,
  seed = 1234
)
```

## Arguments

- per_sample_df:

  Tidy data.frame, typically the output of
  \[nhood_enrichment_per_sample()\], \[cooccur_ratio_per_sample()\],
  \[cooccur_local_per_sample()\], or \[interaction_spot_per_sample()\].

- value:

  Name of the column to test (e.g. "log2_oe", "zscore", "ratio", "mean",
  "spots_per_1k_cells").

- group_key:

  Name of the group column. Defaults to "group".

- patient_key:

  Optional patient column. Required for \`method = "lmm"\` to use as a
  random effect; used as the permutation block for \`method = "perm"\`.
  If \`NULL\` and the data has a \`patient\` column, that column is only
  used to detect pseudoreplication.

- method:

  Statistical test: \* "wilcox" — Wilcoxon rank-sum (two groups); exact
  p-values for small samples without ties, normal approximation
  otherwise. \* "t" — Welch's t-test (two groups) \* "lmm" — linear
  mixed model \`value ~ group + covariates + (1\|patient)\` via
  \`lme4::lmer\`. p-values use Satterthwaite degrees of freedom when the
  \`lmerTest\` package is installed, and otherwise a t reference with
  \`n_patients - n_fixed_effects\` degrees of freedom (group is a
  patient-level factor). Falls back to \`lm()\` when every patient
  contributes a single row. Requires the \`lme4\` package. \* "perm" —
  group-label permutation test on the mean difference (blocked by
  patient if \`patient_key\` is supplied). All relabelings are
  enumerated when there are at most \`n_perms\` of them (exact test).

- n_perms:

  Number of permutations for \`method = "perm"\`.

- adjust:

  Multiple-testing adjustment method passed to \[stats::p.adjust()\].

- pair_keys:

  Column names that together identify a cluster pair. Defaults to
  \`c("cluster_i", "cluster_j")\`. Set to a single column name for cases
  like interaction_spot_per_sample (e.g. "target_cluster").

- ref_group:

  Name of the reference group (e.g. "control"). The effect is reported
  as \`mean_test - mean_ref\`. When \`NULL\`, groups are sorted
  alphabetically and the first is used as reference (a message says
  which); set it explicitly, because e.g. "case" sorts before "control".

- covariates:

  Optional character vector of additional columns (e.g. age, sex, batch,
  \`n_cells\`) included as fixed effects. Only used with \`method =
  "lmm"\`.

- symmetric:

  If \`TRUE\`, test each unordered pair once (rows with \`cluster_i \<=
  cluster_j\`). Neighborhood enrichment scores are (nearly) symmetric,
  so testing both (i, j) and (j, i) doubles the multiple-testing burden
  without adding information.

- min_n_per_group:

  Minimum number of finite observations required in each group; pairs
  below this are skipped.

- seed:

  Random seed for the permutation test.

## Value

A data.frame with the cluster pair columns, group sizes and means,
\`effect\` (test group mean minus reference group mean), for \`method =
"lmm"\` also \`estimate\` (the adjusted model coefficient),
\`statistic\`, raw \`p\`, and \`padj\`. Attributes \`method\`,
\`groups\`, \`value\` and \`p_method\` describe the test.

## Unit of analysis

In a case-control study the independent unit is the \*patient\*, not the
image. With several images per patient either aggregate first (\`unit =
"patient"\` in the \`\*\_per_sample()\` call) and use \`"wilcox"\` /
\`"t"\`, or keep images and use \`"lmm"\` or \`"perm"\` with
\`patient_key\`. \`compare_groups()\` warns when \`"wilcox"\` / \`"t"\`
(or \`"perm"\` without \`patient_key\`) are given several rows per
patient.

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 3,
                          group_close_ratio = list(case = 0.8, control = 0.2),
                          n_types = 4, n_cells = 200, max_loc = 250,
                          test_type = "distribute", distance_param = 8,
                          seed = 1)
ps <- nhood_enrichment_per_sample(df, sample_key = "sample_id",
                                  group_key = "group",
                                  cluster_key = "cell_type",
                                  patient_key = "patient",
                                  neighbors.k = 8, n_perms = 20, n_jobs = 1)
cmp <- compare_groups(ps, value = "log2_oe", method = "wilcox",
                      ref_group = "control", symmetric = TRUE)
head(cmp)
#>     cluster_i   cluster_j n_total n_control n_case mean_control  mean_case
#> 1 cell_type_1 cell_type_3       6         3      3   0.08777052 -0.1548648
#> 2 cell_type_1 cell_type_4       6         3      3   0.07356400 -0.4615872
#> 3 cell_type_2 cell_type_4       6         3      3  -0.07560648 -0.4360321
#> 4 cell_type_4 cell_type_4       6         3      3  -0.06441771  0.1496259
#> 5 cell_type_1 cell_type_2       6         3      3   0.25911989  0.1589639
#> 6 cell_type_2 cell_type_3       6         3      3  -0.05225977 -0.2397878
#>       effect statistic   p      padj
#> 1 -0.2426353         9 0.1 0.2500000
#> 2 -0.5351512         9 0.1 0.2500000
#> 3 -0.3604256         9 0.1 0.2500000
#> 4  0.2140436         0 0.1 0.2500000
#> 5 -0.1001560         8 0.2 0.3333333
#> 6 -0.1875280         8 0.2 0.3333333
```
