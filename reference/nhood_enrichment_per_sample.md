# Per-sample neighborhood enrichment

Run \[nhood_enrichment()\] independently for each sample (image) in the
input, returning a tidy long data.frame with one row per \`sample_id x
cluster_i x cluster_j\`. The cluster factor levels are harmonized across
samples so the resulting table is suitable for downstream group
comparison with \[compare_groups()\].

## Usage

``` r
nhood_enrichment_per_sample(
  obj,
  sample_key,
  group_key,
  cluster_key,
  patient_key = NULL,
  unit = c("image", "patient"),
  cluster_levels = NULL,
  neighbors.k = 30,
  connectivity_key = "nn",
  transformation = FALSE,
  n_perms = 100,
  seed = 1938493,
  n_jobs = 1
)
```

## Arguments

- obj:

  A Seurat object, a list of Seurat objects, or a data.frame with x, y,
  cluster_key, and sample_key columns.

- sample_key:

  Column / image identifier defining a sample.

- group_key:

  Column carrying the disease / condition label.

- cluster_key:

  Column with cluster labels.

- patient_key:

  Optional column with patient ID (random effect / permutation block).

- unit:

  "image" returns one row per image; "patient" averages across images
  within a patient (requires \`patient_key\`).

- cluster_levels:

  Optional character vector of cluster levels to use as common dimnames
  across samples. Defaults to the union across all samples.

- neighbors.k:

  Number of neighbors to use.

- connectivity_key:

  Type of graph: "nn" or "snn".

- transformation:

  If \`TRUE\`, each link from cell u is weighted 1 / (1 + d_u), with d_u
  the number of cells that chose u as a neighbour (a mild down-weighting
  of hub cells). \`FALSE\` (default since 0.99.3) counts every kNN link
  once, as squidpy does. With a kNN graph every cell sends k links at
  any density, so the weighting is not a density correction (density is
  handled by the label shuffles); in simulations it kept the calibration
  but lowered power.

- n_perms:

  Number of permutations.

- seed:

  Random seed.

- n_jobs:

  Number of parallel jobs. \`1\` runs sequentially.

## Value

A data.frame (also tagged with class \`cohaluSample\`) with columns
\`sample_id\`, \`cluster_i\`, \`cluster_j\`, \`zscore\`, \`count\`,
\`expected\`, \`log2_oe\`, \`group\`, \`patient\`, plus \`n_cells\`
(cells in the sample) and \`n_i\` / \`n_j\` (cells of \`cluster_i\` /
\`cluster_j\`).

## Choosing the value to compare

The permutation z-score measures \*statistical evidence\* within one
sample and grows roughly with the square root of the number of cells, so
samples with more cells (larger images, denser tissue) get larger \|z\|
for the same spatial pattern. For between-group comparison, \`log2_oe\`
(log2 observed / permutation-expected) is an effect size that does not
scale with sample size and is usually the better choice. \`count\` is
additionally confounded by cell-type composition and should not be
compared directly.

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
head(ps)
#>   sample_id   cluster_i   cluster_j      zscore count expected      log2_oe
#> 1    case_1 cell_type_1 cell_type_1 -0.05996466    97    97.50 -0.002634965
#> 2    case_1 cell_type_2 cell_type_1  2.51137150   116    91.10  0.353525585
#> 3    case_1 cell_type_3 cell_type_1  1.11862222    90    79.65  0.183330821
#> 4    case_1 cell_type_4 cell_type_1 -0.48086276    55    58.90 -0.084194562
#> 5    case_1 cell_type_1 cell_type_2  2.68011895   112    91.05  0.300985044
#> 6    case_1 cell_type_2 cell_type_2  0.52102225   207   200.40  0.049241298
#>   group patient n_cells n_i n_j
#> 1  case  case_1     200  41  41
#> 2  case  case_1     200  63  41
#> 3  case  case_1     200  55  41
#> 4  case  case_1     200  41  41
#> 5  case  case_1     200  41  63
#> 6  case  case_1     200  63  63
```
