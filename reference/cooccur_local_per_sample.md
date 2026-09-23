# Per-sample local co-occurrence summary

Run \[cooccur_local()\] for one (cluster_x, cluster_y) pair on each
sample and summarize the per-cell scores at the sample level. Output is
a tidy data.frame suitable for group comparison.

## Usage

``` r
cooccur_local_per_sample(
  obj,
  sample_key,
  group_key,
  cluster_key,
  cluster_x,
  cluster_y,
  patient_key = NULL,
  unit = c("image", "patient"),
  neighbors.k = 20,
  radius = 30,
  maxnsteps = 1,
  summarize = c("mean", "q90", "pos_rate")
)
```

## Arguments

- obj:

  A Seurat object, a list of Seurat objects, or a data.frame.

- sample_key, group_key, cluster_key, patient_key:

  Sample / group / cluster / patient column names. See
  \[build_sample_design()\].

- cluster_x:

  First cluster of interest.

- cluster_y:

  Second cluster of interest.

- unit:

  "image" or "patient" (averaged across images of the same patient).

- neighbors.k:

  Number of neighbors.

- radius:

  Radius for neighborhood.

- maxnsteps:

  Maximum number of diffusion steps. Each step computes \`s \<- (A + I)
  D^-1 s\` from the previous step; diffusion stops early once the
  kurtosis of the scores decreases by less than 3 between steps (checked
  after step 3). \`0\` returns the raw 0/1 indicator. Versions \<=
  0.99.1 always performed a single step regardless of \`maxnsteps\`.

- summarize:

  Character vector of summary statistics to compute: "mean", "q90" (90th
  percentile), and / or "pos_rate" (fraction of cells with score \> 0).

## Value

A data.frame with one row per sample carrying the requested summary
statistics, plus \`n_cells\`, \`n_i\` and \`n_j\` (cells of
\`cluster_x\` / \`cluster_y\`).

## Caution

The local score has no permutation null, so its sample-level summaries
increase with the abundance of \`cluster_x\` and \`cluster_y\`. When the
two groups differ in cell-type composition, check \`n_i\` / \`n_j\` (or
adjust for them with \`covariates\` in \[compare_groups()\]) before
interpreting a group difference as a change in co-localization.

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 3,
                          group_close_ratio = list(case = 0.8, control = 0.2),
                          n_types = 4, n_cells = 200, max_loc = 250,
                          test_type = "distribute", distance_param = 8,
                          seed = 1)
cooccur_local_per_sample(df, sample_key = "sample_id", group_key = "group",
                         cluster_key = "cell_type",
                         cluster_x = "cell_type_1", cluster_y = "cell_type_2",
                         patient_key = "patient", neighbors.k = 10,
                         radius = 20)
#>   sample_id   cluster_i   cluster_j  mean       q90 pos_rate n_cells n_i n_j
#> 1    case_1 cell_type_1 cell_type_2 0.530 0.8663170    0.920     200  41  63
#> 2    case_2 cell_type_1 cell_type_2 0.545 0.8285827    0.960     200  49  52
#> 3    case_3 cell_type_1 cell_type_2 0.525 0.8651515    0.950     200  51  44
#> 4 control_1 cell_type_1 cell_type_2 0.335 0.6922113    0.870     200  42  42
#> 5 control_2 cell_type_1 cell_type_2 0.395 0.6866178    0.935     200  53  52
#> 6 control_3 cell_type_1 cell_type_2 0.270 0.6303030    0.750     200  53  42
#>     group   patient
#> 1    case    case_1
#> 2    case    case_2
#> 3    case    case_3
#> 4 control control_1
#> 5 control control_2
#> 6 control control_3
```
