test_that("generate_sim_groups produces a multi-sample data.frame", {
  df <- generate_sim_groups(
    n_samples_per_group = 2,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 5, n_cells = 200, max_loc = 400,
    distance_param = 20, between_sample_noise = 0,
    seed = 42
  )
  expect_s3_class(df, "data.frame")
  expect_true(all(c("x", "y", "cell_type", "sample_id", "group", "patient") %in% colnames(df)))
  expect_setequal(unique(df$group), c("disease", "control"))
  expect_length(unique(df$sample_id), 4L)
})

test_that("build_sample_design extracts sample/group from a data.frame", {
  df <- generate_sim_groups(
    n_samples_per_group = 2,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 4, n_cells = 100, max_loc = 300,
    distance_param = 20, between_sample_noise = 0, seed = 1
  )
  ds <- build_sample_design(df, sample_key = "sample_id", group_key = "group",
                            patient_key = "patient")
  expect_s3_class(ds, "data.frame")
  expect_length(unique(ds$sample_id), 4L)
  expect_setequal(unique(ds$group), c("disease", "control"))
  expect_true(all(!is.na(ds$patient)))
})

test_that("nhood_enrichment_per_sample returns tidy per-sample results", {
  skip_on_cran()
  df <- generate_sim_groups(
    n_samples_per_group = 2,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 4, n_cells = 200, max_loc = 400,
    distance_param = 20, between_sample_noise = 0, seed = 7
  )
  res <- nhood_enrichment_per_sample(
    df, sample_key = "sample_id", group_key = "group",
    cluster_key = "cell_type", patient_key = "patient",
    neighbors.k = 10, n_perms = 20, n_jobs = 1, seed = 7
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("sample_id", "cluster_i", "cluster_j",
                    "zscore", "count", "group", "patient") %in% colnames(res)))
  expect_length(unique(res$sample_id), 4L)
  # 4 cluster levels -> 16 pairs per sample
  expect_equal(nrow(res), 4L * 4L * 4L)
})

test_that("compare_groups detects the planted disease effect", {
  skip_on_cran()
  # Use test_type = "distribute" because it places type-1 / type-2 cells at
  # pairwise distance ~distance_param, which is the regime the local
  # neighborhood graph (k = 10) is sensitive to. The "circle" simulation
  # produces co-occurrence at a coarser scale (inner disk vs outer ring)
  # and tends to look like local segregation under small-k neighborhoods.
  df <- generate_sim_groups(
    n_samples_per_group = 3,
    group_close_ratio = list(disease = 0.9, control = 0.1),
    n_types = 4, n_cells = 400, max_loc = 400,
    test_type = "distribute",
    distance_param = 15, between_sample_noise = 0.02, seed = 11
  )
  res <- nhood_enrichment_per_sample(
    df, sample_key = "sample_id", group_key = "group",
    cluster_key = "cell_type", patient_key = "patient",
    neighbors.k = 10, n_perms = 30, n_jobs = 1, seed = 11
  )
  cmp <- compare_groups(res, value = "zscore", group_key = "group",
                        method = "wilcox", adjust = "BH",
                        ref_group = "control")
  expect_s3_class(cmp, "data.frame")
  expect_true(all(c("cluster_i", "cluster_j", "effect", "p", "padj") %in% colnames(cmp)))

  # Cluster pair 1-2 was planted as strongly co-localized in disease.
  # With ref_group = "control", effect = mean_disease - mean_control > 0.
  target <- cmp[cmp$cluster_i == "cell_type_1" & cmp$cluster_j == "cell_type_2", ]
  expect_gt(nrow(target), 0L)
  expect_gt(target$effect[1], 0)
})

test_that("cooccur_ratio_per_sample returns ratio + count per pair", {
  skip_on_cran()
  df <- generate_sim_groups(
    n_samples_per_group = 2,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 4, n_cells = 200, max_loc = 400,
    distance_param = 20, between_sample_noise = 0, seed = 21
  )
  res <- cooccur_ratio_per_sample(
    df, sample_key = "sample_id", group_key = "group",
    cluster_key = "cell_type", patient_key = "patient",
    radius = 30, k = 10
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("ratio", "count") %in% colnames(res)))
  expect_length(unique(res$sample_id), 4L)
})

test_that("cooccur_local_per_sample summarizes per sample", {
  skip_on_cran()
  df <- generate_sim_groups(
    n_samples_per_group = 2,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 4, n_cells = 200, max_loc = 400,
    distance_param = 20, between_sample_noise = 0, seed = 31
  )
  res <- cooccur_local_per_sample(
    df, sample_key = "sample_id", group_key = "group",
    cluster_key = "cell_type",
    cluster_x = "cell_type_1", cluster_y = "cell_type_2",
    patient_key = "patient",
    neighbors.k = 10, radius = 30, maxnsteps = 1
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean", "q90", "pos_rate") %in% colnames(res)))
  expect_length(unique(res$sample_id), 4L)
})

test_that("compare_groups with method = 'perm' returns finite p-values", {
  skip_on_cran()
  df <- generate_sim_groups(
    n_samples_per_group = 3,
    group_close_ratio = list(disease = 0.8, control = 0.2),
    n_types = 4, n_cells = 200, max_loc = 400,
    distance_param = 20, between_sample_noise = 0.02, seed = 41
  )
  res <- cooccur_ratio_per_sample(
    df, sample_key = "sample_id", group_key = "group",
    cluster_key = "cell_type", patient_key = "patient",
    radius = 30, k = 10
  )
  cmp <- compare_groups(res, value = "ratio", group_key = "group",
                        patient_key = "patient",
                        method = "perm", n_perms = 100, seed = 41)
  expect_s3_class(cmp, "data.frame")
  expect_true(all(cmp$p >= 0 & cmp$p <= 1))
})
