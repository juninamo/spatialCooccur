# Convert simulated cells to a Seurat object with one FOV per sample

Build a minimal Seurat object from a coordinate table such as the output
of \[generate_sim()\] or \[generate_sim_groups()\]: every sample becomes
a centroid-based FOV named after it, and all other columns are stored in
\`meta.data\`. The expression matrix is a small placeholder, so the
object is meant for exercising the Seurat input path of the spatial
functions, not for expression analysis.

## Usage

``` r
sim_to_seurat(df, sample_key = "sample_id", n_features = 5)
```

## Arguments

- df:

  A data.frame with \`x\`, \`y\` and one row per cell.

- sample_key:

  Column identifying the sample (image). If absent, all cells are put in
  a single FOV named \`"fov"\`.

- n_features:

  Number of placeholder features in the count matrix.

## Value

A Seurat object with a \`cell\` column in \`meta.data\` and one image
per sample in \`@images\`.

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 2, n_types = 4,
                          n_cells = 150, max_loc = 250,
                          test_type = "distribute", distance_param = 10)
seu <- sim_to_seurat(df)
SeuratObject::Images(seu)
#> [1] "disease_1" "disease_2" "control_1" "control_2"
```
