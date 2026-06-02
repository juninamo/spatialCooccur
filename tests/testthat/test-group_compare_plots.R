test_that("plot helpers return ggplot objects", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  df <- generate_sim_groups(
    n_samples_per_group = 2,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 4, n_cells = 200, max_loc = 400,
    test_type = "distribute",
    distance_param = 15, between_sample_noise = 0.02, seed = 99
  )
  per_sample <- nhood_enrichment_per_sample(
    df, sample_key = "sample_id", group_key = "group",
    cluster_key = "cell_type", patient_key = "patient",
    neighbors.k = 10, n_perms = 30, n_jobs = 1, seed = 99
  )
  cmp <- compare_groups(per_sample, value = "zscore",
                        method = "wilcox", adjust = "BH",
                        ref_group = "control")

  expect_s3_class(plot_group_delta_heatmap(cmp), "ggplot")
  expect_s3_class(plot_volcano_groups(cmp, label_top = 3), "ggplot")
  expect_s3_class(
    plot_pair_boxplot(per_sample, value = "zscore",
                      pairs = data.frame(cluster_i = "cell_type_1",
                                         cluster_j = "cell_type_2"),
                      add_p = TRUE,
                      ref_group = "control"),
    "ggplot"
  )
})
