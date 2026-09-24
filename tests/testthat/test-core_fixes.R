test_that("compute_count matches a brute-force submatrix sum", {
  set.seed(1)
  n <- 60
  adj <- Matrix::rsparsematrix(n, n, density = 0.1, rand.x = function(k) runif(k))
  cl <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  lab <- paste0("Cluster", cl)
  got <- compute_count(adj, lab, lab, 3, cl, transformation = TRUE)
  for (i in levels(cl)) for (j in levels(cl)) {
    expect_equal(got[paste0("Cluster", i), paste0("Cluster", j)],
                 sum(adj[cl == i, cl == j]))
  }
})

test_that("nhood_enrichment is reproducible for a given seed", {
  df <- generate_sim(close_ratio = 0.5, n_types = 4, n_cells = 300,
                     max_loc = 300, test_type = "distribute",
                     distance_param = 10, seed = 3)
  a <- nhood_enrichment(df, "cell_type", neighbors.k = 8, n_perms = 30, seed = 9, n_jobs = 1)
  b <- nhood_enrichment(df, "cell_type", neighbors.k = 8, n_perms = 30, seed = 9, n_jobs = 1)
  expect_identical(a, b)
  expect_true(all(c("zscore", "count", "expected", "log2_oe") %in% names(a)))
})

test_that("nhood_enrichment z-scores are calibrated under spatial randomness", {
  skip_on_cran()
  set.seed(10)
  z <- replicate(4, {
    n <- 1200
    d <- data.frame(x = runif(n, 0, 500), y = runif(n, 0, 500),
                    cell_type = factor(sample(paste0("t", 1:5), n, replace = TRUE)))
    rownames(d) <- paste0("c", seq_len(n))
    r <- nhood_enrichment(d, "cell_type", neighbors.k = 10, n_perms = 100,
                          seed = sample.int(1e6, 1), n_jobs = 1)$zscore
    c(diag = mean(diag(r)), off = mean(r[upper.tri(r)]))
  })
  # Independent row / column shuffling gave diag ~ +11 and off-diag ~ -2.
  expect_lt(abs(mean(z["diag", ])), 1)
  expect_lt(abs(mean(z["off", ])), 1)
})

test_that("absent cell types give NA rather than dropping rows", {
  df <- generate_sim(close_ratio = 0.5, n_types = 4, n_cells = 200,
                     max_loc = 300, test_type = "distribute",
                     distance_param = 10, seed = 3)
  df$cell_type <- factor(as.character(df$cell_type),
                         levels = c(levels(df$cell_type), "absent"))
  r <- nhood_enrichment(df, "cell_type", neighbors.k = 8, n_perms = 10, seed = 1, n_jobs = 1)
  expect_equal(dim(r$zscore), c(5L, 5L))
  expect_true(all(is.na(r$zscore["Clusterabsent", ])))
})

test_that("seeded functions leave the caller's RNG stream untouched", {
  d <- data.frame(sample_id = paste0("s", 1:8), patient = paste0("s", 1:8),
                  group = rep(c("case", "control"), each = 4),
                  cluster_i = "A", cluster_j = "B", value = rnorm(8))
  set.seed(42); a1 <- runif(1); a2 <- runif(1)
  set.seed(42); b1 <- runif(1)
  invisible(compare_groups(d, value = "value", method = "perm", ref_group = "control", n_perms = 50))
  b2 <- runif(1)
  expect_identical(c(a1, a2), c(b1, b2))
})

test_that("cooccur_local diffusion iterates from the previous step", {
  df <- generate_sim(close_ratio = 0.8, n_types = 5, n_cells = 400, max_loc = 300,
                     test_type = "distribute", distance_param = 10, seed = 2)
  s0 <- cooccur_local(df, "cell_type_1", "cell_type_2", neighbors.k = 10, radius = 20, maxnsteps = 0)[[1]]
  s1 <- cooccur_local(df, "cell_type_1", "cell_type_2", neighbors.k = 10, radius = 20, maxnsteps = 1)[[1]]
  s3 <- cooccur_local(df, "cell_type_1", "cell_type_2", neighbors.k = 10, radius = 20, maxnsteps = 3)[[1]]
  expect_true(all(s0 %in% c(0, 1)))
  expect_false(isTRUE(all.equal(s1, s3)))      # extra steps now change the result
  expect_equal(sum(s3), sum(s0))               # diffusion conserves total mass
  expect_lt(sd(s3), sd(s1))                    # and keeps smoothing
})

test_that("generate_sim keeps relocated cells inside the tissue", {
  for (dp in c(10, 100)) {
    d <- generate_sim(close_ratio = 0.5, n_types = 8, n_cells = 1500, max_loc = 600,
                      test_type = "distribute", distance_param = dp, seed = 3)
    expect_true(all(d$x >= 0 & d$x <= 600 & d$y >= 0 & d$y <= 600))
  }
})

test_that("centred log2_oe is unbiased for rare cell types under the null", {
  set.seed(11)
  est <- sapply(1:6, function(i) {
    n <- 1200
    d <- data.frame(x = runif(n, 0, 550), y = runif(n, 0, 550),
                    cell_type = sample(paste0("t", 1:20), n, TRUE))
    rownames(d) <- paste0("c", seq_len(n))
    r <- nhood_enrichment(d, "cell_type", neighbors.k = 10, n_perms = 60, seed = i, n_jobs = 1)
    up <- upper.tri(r$log2_oe)
    c(raw = mean(r$log2_oe_raw[up]), centred = mean(r$log2_oe[up]))
  })
  expect_lt(mean(est["raw", ]), -0.02)            # the uncorrected log ratio is biased
  expect_lt(abs(mean(est["centred", ])), 0.01)    # the centred one is not
})

test_that("pvalue / padj are symmetric and max-T adjusted", {
  d <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 1000, max_loc = 500,
                    test_type = "distribute", distance_param = 10, seed = 1)
  rownames(d) <- paste0("c", seq_len(nrow(d)))
  r <- nhood_enrichment(d, "cell_type", neighbors.k = 10, n_perms = 99, seed = 1, n_jobs = 1)
  expect_true(isSymmetric(unname(r$padj)))
  expect_true(all(r$padj >= r$pvalue - 1e-12 | r$padj >= 1 / 100, na.rm = TRUE))
  expect_gte(min(r$padj, na.rm = TRUE), 1 / 100)
  expect_lt(r$padj[1, 2], 0.05)                    # the planted pair
})

test_that("directional contact / dominance statistics separate the two points of view", {
  set.seed(1); n <- 2000
  d <- data.frame(x = runif(n, 0, 500), y = runif(n, 0, 500), cell_type = sample(c("B", "O"), n, TRUE, prob = c(0.4, 0.6)))
  ctr <- cbind(runif(60, 20, 480), runif(60, 20, 480))
  d <- rbind(d, data.frame(x = ctr[, 1], y = ctr[, 2], cell_type = "A"),
             data.frame(x = rep(ctr[, 1], each = 6) + rnorm(360, 0, 6), y = rep(ctr[, 2], each = 6) + rnorm(360, 0, 6), cell_type = "B"))
  rownames(d) <- paste0("c", seq_len(nrow(d))); d$cell_type <- factor(d$cell_type)
  r <- nhood_enrichment(d, cluster_key = "cell_type", neighbors.k = 10, n_perms = 99, seed = 1, n_jobs = 1)
  nm <- function(m) { dimnames(m) <- lapply(dimnames(m), function(v) sub("^Cluster", "", v)); m }
  dom <- nm(r$dominance_log2_oe); con <- nm(r$contact_log2_oe)
  expect_gt(dom["A", "B"], 0.5)                # A cells are surrounded by B
  expect_true(is.na(dom["B", "A"]) || abs(dom["B", "A"]) < 0.2)   # B neighbourhoods are not dominated by A
  expect_gt(con["B", "A"], 0.3)                # more B cells than expected touch an A cell
  expect_false(isSymmetric(unname(round(nm(r$dominance), 3))))
  expect_true(all(r$contact_padj >= 1 / 100, na.rm = TRUE))
})
