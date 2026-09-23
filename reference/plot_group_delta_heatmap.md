# Heatmap of the disease-group effect for every cluster pair

Plot \`effect\` (or another column from a \[compare_groups()\] result)
as a cluster_i x cluster_j heatmap, with optional significance markers.

## Usage

``` r
plot_group_delta_heatmap(
  compare_df,
  value = "effect",
  significance = "padj",
  sig_threshold = 0.05,
  pair_keys = c("cluster_i", "cluster_j"),
  palette = c("RdBu", "viridis"),
  limits = NULL
)
```

## Arguments

- compare_df:

  Output of \[compare_groups()\].

- value:

  Column to map to fill. Defaults to "effect".

- significance:

  Column name carrying p-values to overlay as labels (set to NULL to
  skip).

- sig_threshold:

  Adjusted p-value threshold for the asterisk marker.

- pair_keys:

  Two column names identifying the cluster pair on the x and y axes.
  Defaults to \`c("cluster_i", "cluster_j")\`.

- palette:

  One of "RdBu" (diverging, default) or "viridis".

- limits:

  Optional length-2 numeric vector for the fill scale limits. If NULL,
  symmetric limits around 0 are used for diverging palettes.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
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
  plot_group_delta_heatmap(cmp)
}
```
