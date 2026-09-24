test_that(".grf_fft gives unit-variance fields with an RBF covariance", {
  set.seed(1)
  f <- spatialCooccur:::.grf_fft(200, 200, delta = 1, ell = 5)
  expect_equal(var(as.vector(f)), 1, tolerance = 0.15)
  lag_cor <- cor(as.vector(f[1:195, ]), as.vector(f[6:200, ]))
  expect_equal(lag_cor, exp(-0.5), tolerance = 0.08)
})

test_that("bin_transcripts keeps every transcript", {
  tx <- simulate_transcripts(size = 60, rate = 0.02, n_genes_per_set = 2, seed = 2)
  b <- bin_transcripts(tx, bin_size = 5)
  expect_s3_class(b, "binned_transcripts")
  expect_equal(sum(b$counts), nrow(tx))
  expect_equal(nrow(b$counts), b$grid$nx * b$grid$ny)
  expect_setequal(b$genes, unique(tx$gene))
})

test_that("pcf_cross is ~0 on the log scale under complete spatial randomness", {
  set.seed(3)
  n <- 6000
  tx <- data.frame(x = runif(n, 0, 300), y = runif(n, 0, 300),
                   gene = sample(c("a", "b"), n, replace = TRUE))
  b <- bin_transcripts(tx, bin_size = 5, tissue_radius = Inf)
  ab <- pcf_cross(b, "a", "b", r_max = 40)
  aa <- pcf_cross(b, "a", r_max = 40)          # self-pairs must be removed
  expect_lt(max(abs(ab$log_g)), 0.15)
  expect_lt(max(abs(aa$log_g)), 0.15)
})

test_that("pcf_cross and the model recover planted co-localization", {
  skip_on_cran()
  tx <- simulate_transcripts(size = 200, rate = 0.03, n_genes_per_set = 3,
                             coloc = c(A = 1.2, B = 1.2), set_sd = 0.5, seed = 4)
  sets <- with(unique(tx[, c("gene", "gene_set")]), split(gene, gene_set))
  b <- bin_transcripts(tx, bin_size = 5, tissue_radius = Inf)
  ab <- pcf_cross(b, sets$A, sets$B, r_max = 20)
  ac <- pcf_cross(b, sets$A, sets$C, r_max = 20)
  expect_gt(ab$log_g[2], ac$log_g[2] + 0.3)
  fit <- fit_spatial_rff(b, n_factors = 3, n_features = 32, max_iter = 150)
  expect_s3_class(fit, "spatial_rff_fit")
  expect_equal(dim(fit$L), c(9L, 3L))
  m_ab <- rff_pair_correlation(fit, sets$A, sets$B, r = 10, type = "composition")$log_g
  m_ac <- rff_pair_correlation(fit, sets$A, sets$C, r = 10, type = "composition")$log_g
  expect_gt(m_ab, m_ac)
  expect_gt(m_ab, 0.15)
  g_ab <- rff_pair_correlation(fit, sets$A, sets$B, r = 10, type = "composition",
                               method = "gaussian")$log_g
  expect_gt(g_ab, 0.3)
})

test_that("colocalization_per_sample feeds compare_groups", {
  skip_on_cran()
  tx <- simulate_transcripts_groups(n_samples_per_group = 3, size = 100, rate = 0.03,
                                    n_genes_per_set = 2, seed = 5)
  sets <- with(unique(tx[, c("gene", "gene_set")]), split(gene, gene_set))
  res <- colocalization_per_sample(tx, "sample_id", "group", "patient", sets,
                                   pairs = list(c("A", "B")), r = c(8, 16))
  expect_equal(nrow(res), 6L * 2L)
  expect_true(all(c("cluster_i", "cluster_j", "r", "log_g", "group") %in% colnames(res)))
  cmp <- compare_groups(res, value = "log_g", method = "wilcox", ref_group = "control",
                        pair_keys = c("cluster_i", "cluster_j", "r"))
  expect_equal(nrow(cmp), 2L)
})

test_that("fit_spatial_rff supports a smoothed-total offset", {
  tx <- simulate_transcripts(size = 80, rate = 0.03, n_genes_per_set = 2, seed = 6)
  b <- bin_transcripts(tx, bin_size = 5, tissue_radius = Inf)
  fit <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 30,
                         offset = "smoothed_total", offset_genes = c("C_1", "C_2"))
  expect_false(fit$has_density)
  expect_identical(fit$offset, "smoothed_total")
  expect_error(fit_spatial_rff(b, offset = "smoothed_total", offset_genes = "nope"),
               "offset_genes")
  lg <- rff_pair_correlation(fit, c("A_1", "A_2"), c("B_1", "B_2"), r = 10, type = "composition")
  expect_true(is.finite(lg$log_g))
})

test_that("pcf_matrix agrees with pcf_cross", {
  tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 2, seed = 12)
  sets <- split(attr(tx, "truth")$genes, attr(tx, "truth")$set_of)
  b <- bin_transcripts(tx, bin_size = 5)
  m <- pcf_matrix(b, sets, r_max = 30)
  expect_equal(nrow(m), 6L * length(unique(m$r)))
  ab <- pcf_cross(b, sets$A, sets$B, r_max = 30)
  ab_rel <- pcf_cross(b, sets$A, sets$B, r_max = 30, relative = TRUE)
  aa <- pcf_cross(b, sets$A, r_max = 30)
  mm <- m[m$cluster_i == "A" & m$cluster_j == "B", ]
  expect_equal(mm$log_g, ab$log_g, tolerance = 1e-8)
  expect_equal(mm$log_g_rel, ab_rel$log_g, tolerance = 1e-8)
  expect_equal(m$log_g[m$cluster_i == "A" & m$cluster_j == "A"], aa$log_g, tolerance = 1e-8)
})

test_that("pcf_matrix accepts any labelled points, e.g. pixel-level factors", {
  set.seed(1); n <- 3000
  cx <- runif(15, 0, 300); cy <- runif(15, 0, 300); k <- sample(15, n, TRUE)
  px <- data.frame(X = c(cx[k] + rnorm(n, 0, 8), runif(n, 0, 300)),
                   Y = c(cy[k] + rnorm(n, 0, 8), runif(n, 0, 300)),
                   K1 = c(sample(c("1", "2"), n, TRUE), rep("3", n)))
  b <- bin_transcripts(px, bin_size = 4, x_col = "X", y_col = "Y", gene_col = "K1")
  pm <- pcf_matrix(b, list(F1 = "1", F2 = "2", F3 = "3"), r_max = 30)
  short <- pm[pm$r > 0 & pm$r <= 20 & is.finite(pm$log_g_rel), ]
  expect_gt(mean(short$log_g_rel[short$cluster_i == "F1" & short$cluster_j == "F2"]), 0.5)
  expect_lt(mean(short$log_g_rel[short$cluster_i == "F1" & short$cluster_j == "F3"]), -0.3)
})
