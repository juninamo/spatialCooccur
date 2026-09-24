test_that("rff_fields() reproduces the fitted grid and averages per cell", {
  tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3, seed = 2)
  b <- bin_transcripts(tx, bin_size = 6)
  fit <- fit_spatial_rff(b, n_factors = 2, n_features = 24, max_iter = 30)
  keep <- fit$in_tissue
  ctr <- b$coords[keep, c("x", "y")]
  f <- rff_fields(fit, ctr)
  expect_equal(unname(as.matrix(f)), unname(fit$field_grid[keep, , drop = FALSE]), tolerance = 1e-8)
  expect_equal(colnames(f), c("density", "factor1", "factor2"))
  # fits made before `center` was stored give the same values
  old <- fit; old$center <- NULL
  expect_equal(rff_fields(old, ctr), f, tolerance = 1e-10)
  # per-cell averages
  tx$cell_id <- paste(floor(tx$x / 20), floor(tx$y / 20))
  pc <- rff_fields(fit, tx, by = "cell_id", which = "factors")
  pt <- rff_fields(fit, tx, which = "factors")
  one <- pc$cell_id[1]
  expect_equal(pc$factor1[1], mean(pt$factor1[tx$cell_id == one]), tolerance = 1e-10)
  expect_equal(pc$n_points[1], sum(tx$cell_id == one))
  expect_equal(sum(pc$n_points), nrow(tx))
  # points far outside the tissue are flagged
  expect_warning(rff_fields(fit, data.frame(x = 1e4, y = 1e4)), "outside the fitted tissue")
  expect_error(rff_fields(fit, data.frame(a = 1)), "Columns not found")
})
