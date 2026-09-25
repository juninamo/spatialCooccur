test_that("a zero matrix offset reproduces the area-only fit", {
  tx <- simulate_transcripts(size = 80, rate = 0.02, n_genes_per_set = 2, seed = 3)
  b <- bin_transcripts(tx, bin_size = 8)
  f1 <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 20)
  z <- matrix(0, sum(b$coords$in_tissue), length(b$genes), dimnames = list(NULL, b$genes))
  f2 <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 20, offset = z)
  expect_equal(f1$objective, f2$objective, tolerance = 1e-6)
  expect_identical(f2$offset, "matrix")
  expect_error(fit_spatial_rff(b, offset = z[-1, ]), "one row per bin")
})

test_that("rff_offset() absorbs known structure and ARD shrinks unneeded factors", {
  set.seed(4); n <- 30; res <- 6
  dom <- cohalu:::.grf_fft(n, n, res, 30) > 0
  J <- 12; prof <- matrix(rnorm(J * 2, 0, 1.2), J, 2)
  tx <- do.call(rbind, lapply(seq_len(J), function(j) {
    cnt <- rpois(n * n, exp(log(0.1) + ifelse(dom, prof[j, 1], prof[j, 2])) * res^2); idx <- rep(seq_len(n * n), cnt)
    data.frame(x = ((idx - 1) %% n + runif(length(idx))) * res, y = ((idx - 1) %/% n + runif(length(idx))) * res, gene = paste0("g", j))
  }))
  b <- bin_transcripts(tx, bin_size = res, tissue_radius = Inf)
  off <- rff_offset(b, data.frame(domain = factor(as.vector(dom))))
  expect_equal(dim(off), c(n * n, J))
  # the offset carries the domain contrast of each gene
  d1 <- colMeans(off[as.vector(dom), ]) - colMeans(off[!as.vector(dom), ])
  expect_gt(cor(d1, (prof[, 1] - prof[, 2])[match(colnames(off), paste0("g", seq_len(J)))]), 0.95)
  plain <- fit_spatial_rff(b, n_factors = 3, n_features = 32, max_iter = 60)
  resid <- fit_spatial_rff(b, n_factors = 3, n_features = 32, max_iter = 60, offset = off, ard = 30)
  expect_lt(max(resid$factor_strength), 0.5 * max(plain$factor_strength))   # nothing left to explain
  expect_length(resid$factor_strength, 3)
  expect_true(!is.null(resid$settings))
})

test_that("rff_factor_test() returns a p-value per factor", {
  tx <- simulate_transcripts(size = 70, rate = 0.02, n_genes_per_set = 2, seed = 5)
  b <- bin_transcripts(tx, bin_size = 8)
  f <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 20, ard = 2)
  tt <- rff_factor_test(f, b, n_boot = 2, max_iter = 10)
  expect_equal(nrow(tt), 2)
  expect_true(all(tt$p > 0 & tt$p <= 1))
  expect_length(attr(tt, "null_max"), 2)
  expect_true(all(c("program_strength", "uniform_share") %in% names(tt)))
  expect_true(all(tt$program_strength <= tt$strength + 1e-12))
  expect_equal(unname(f$program_strength), unname(sqrt(colSums(sweep(f$L, 2, colMeans(f$L))^2))))
  expect_length(attr(rff_factor_test(f, b, n_boot = 2, max_iter = 10, statistic = "strength"), "null_max"), 2)
})
