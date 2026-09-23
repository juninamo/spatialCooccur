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

# ---- statistical behaviour of compare_groups ------------------------------

.toy_scores <- function(n_pat = 4, n_img = 1, shift = 0, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(c("control", "case"), function(g) {
    do.call(rbind, lapply(seq_len(n_pat), function(p) {
      pe <- rnorm(1, sd = 1)
      data.frame(sample_id = paste0(g, p, "_", seq_len(n_img)),
                 patient = paste0(g, p), group = g,
                 cluster_i = "A", cluster_j = "B",
                 value = pe + rnorm(n_img, sd = 0.2) + (g == "case") * shift)
    }))
  }))
}

test_that("wilcox uses exact p-values for small samples", {
  d <- .toy_scores(n_pat = 3, shift = 10)
  cmp <- compare_groups(d, value = "value", method = "wilcox", ref_group = "control")
  expect_equal(cmp$p, 0.1)
})

test_that("pseudoreplication triggers a warning", {
  d <- .toy_scores(n_pat = 3, n_img = 3)
  expect_warning(compare_groups(d, value = "value", method = "wilcox", ref_group = "control"),
                 "pseudoreplication")
  expect_silent(compare_groups(d, value = "value", method = "perm", patient_key = "patient",
                               ref_group = "control", n_perms = 100))
})

test_that("blocked permutation is exact and respects patient blocks", {
  d <- .toy_scores(n_pat = 3, n_img = 3, shift = 10)
  cmp <- compare_groups(d, value = "value", method = "perm", patient_key = "patient",
                        ref_group = "control", n_perms = 1000)
  # choose(6, 3) = 20 relabelings, two-sided -> smallest p is 2 / 20
  expect_equal(cmp$p, 0.1)
})

test_that("lmm supports covariates and reports an adjusted estimate", {
  skip_if_not_installed("lme4")
  d <- .toy_scores(n_pat = 6, n_img = 3, shift = 2, seed = 4)
  d$batch <- rep(c("b1", "b2"), length.out = nrow(d))
  cmp <- compare_groups(d, value = "value", method = "lmm", patient_key = "patient",
                        covariates = "batch", ref_group = "control")
  expect_true(all(c("estimate", "statistic", "p") %in% colnames(cmp)))
  expect_gt(cmp$estimate, 0)
  expect_true(cmp$p > 0 && cmp$p < 1)
})

test_that("symmetric = TRUE keeps each unordered pair once", {
  d <- .toy_scores(n_pat = 3)
  d2 <- d; d2$cluster_i <- "B"; d2$cluster_j <- "A"
  cmp <- compare_groups(rbind(d, d2), value = "value", method = "t",
                        ref_group = "control", symmetric = TRUE)
  expect_equal(nrow(cmp), 1L)
})

test_that("ref_group = NULL announces the chosen reference", {
  d <- .toy_scores(n_pat = 3)
  expect_message(compare_groups(d, value = "value", method = "t"), "reference")
})

test_that("generate_sim_groups supports several images per patient", {
  df <- generate_sim_groups(n_samples_per_group = 2, n_images_per_patient = 3,
                            group_close_ratio = list(case = 0.5, control = 0.1),
                            n_types = 4, n_cells = 100, max_loc = 200,
                            test_type = "distribute", distance_param = 10, seed = 2)
  expect_length(unique(df$sample_id), 12L)
  expect_length(unique(df$patient), 4L)
})

test_that("nhood_enrichment_per_sample returns log2_oe and composition", {
  skip_on_cran()
  df <- generate_sim_groups(n_samples_per_group = 2,
                            group_close_ratio = list(case = 0.8, control = 0.1),
                            n_types = 4, n_cells = 200, max_loc = 300,
                            test_type = "distribute", distance_param = 8, seed = 3)
  res <- nhood_enrichment_per_sample(df, sample_key = "sample_id", group_key = "group",
                                     cluster_key = "cell_type", patient_key = "patient",
                                     neighbors.k = 8, n_perms = 20, n_jobs = 1)
  expect_true(all(c("log2_oe", "expected", "n_cells", "n_i", "n_j") %in% colnames(res)))
  one <- res[res$sample_id == res$sample_id[1] & res$cluster_i == "cell_type_1", ]
  expect_equal(unique(one$n_i), sum(df$sample_id == res$sample_id[1] & df$cell_type == "cell_type_1"))
})
