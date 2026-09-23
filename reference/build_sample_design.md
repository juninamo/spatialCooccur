# Build a sample design table for disease-group comparisons

Construct a data.frame mapping image / sample identifiers to a disease
group label (and optionally a patient ID). This is the shared metadata
table consumed by the \`\*\_per_sample\` helpers and by
\`compare_groups()\`.

## Usage

``` r
build_sample_design(obj, sample_key, group_key, patient_key = NULL)
```

## Arguments

- obj:

  A Seurat object, a list of Seurat objects, or a data.frame. For Seurat
  input, \`sample_key\` is expected to be a \`meta.data\` column whose
  values match the image names in \`obj@images\`. For data.frame input,
  \`sample_key\` is a column of the data.frame.

- sample_key:

  Name of the column / image identifier.

- group_key:

  Name of the column carrying the disease group (or any condition
  label).

- patient_key:

  Optional column name for patient identifier, used as a random effect
  (or permutation block) downstream. NULL to omit.

## Value

A data.frame with columns \`sample_id\`, \`group\`, \`patient\`, and
\`source_index\` (which list element of the input the sample came from).

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 3,
                          group_close_ratio = list(case = 0.8, control = 0.2),
                          n_types = 4, n_cells = 200, max_loc = 250,
                          test_type = "distribute", distance_param = 8,
                          seed = 1)
build_sample_design(df, sample_key = "sample_id", group_key = "group",
                    patient_key = "patient")
#>           sample_id   group   patient source_index
#> case_1       case_1    case    case_1            1
#> case_2       case_2    case    case_2            1
#> case_3       case_3    case    case_3            1
#> control_1 control_1 control control_1            1
#> control_2 control_2 control control_2            1
#> control_3 control_3 control control_3            1
```
