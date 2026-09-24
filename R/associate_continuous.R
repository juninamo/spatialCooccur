# ---- Association with a continuous sample- or patient-level variable ------
#
# compare_groups() tests a difference between groups of patients. Many
# clinical variables are continuous (CRP, disease activity, number of
# affected joints, age). associate_continuous() tests, for every cell-type
# pair, whether a per-image or per-patient co-localization score changes with
# such a variable, with the patient as the unit of analysis.

#' Associate co-localization with a continuous clinical variable
#'
#' For every cell-type pair (or pair x distance), tests whether a
#' co-localization score from a `*_per_sample()` helper changes with a
#' continuous variable such as CRP, a disease-activity score or age. Images
#' of the same patient are averaged (`"spearman"`, `"lm"`, `"perm"`) or
#' modelled with a random patient intercept (`"lmm"`), so the patient is the
#' unit of analysis.
#'
#' * `"spearman"`: Spearman correlation between patient means and the
#'   variable (exact p-value for small cohorts without ties). Robust to
#'   outliers and to non-linear but monotone relationships.
#' * `"lm"`: linear regression of patient means on the variable, optionally
#'   adjusted for `covariates`; reports the slope.
#' * `"lmm"`: linear mixed model on images,
#'   `score ~ variable + covariates + (1 | patient)` (lmerTest,
#'   Satterthwaite df); uses every image when patients have several.
#' * `"perm"`: Spearman correlation of patient means with a permutation
#'   p-value (the variable is shuffled between patients).
#'
#' @param per_sample_df Output of a `*_per_sample()` helper (one row per
#'   sample x pair), with a `patient` column when patients have several
#'   images.
#' @param variable Name of the numeric column with the clinical variable, or
#'   a named numeric vector indexed by patient (or sample) ID.
#' @param value Score column to test, e.g. `"log2_oe"` or `"log_g_rel"`.
#' @param method One of `"spearman"`, `"lm"`, `"lmm"`, `"perm"`.
#' @param patient_key Column identifying patients (default `"patient"`).
#' @param covariates Optional covariate columns for `"lm"` / `"lmm"`.
#' @param pair_keys Columns identifying one test (default `cluster_i`,
#'   `cluster_j`, plus `r` when present).
#' @param n_perms Permutations for `method = "perm"`.
#' @param adjust Multiple-testing correction passed to [stats::p.adjust()].
#' @param min_n Minimum number of patients with both values.
#' @param seed Random seed for `"perm"`.
#'
#' @return A data.frame with one row per pair: `n_patients`, `estimate`
#'   (Spearman rho, or the slope per unit of `variable`), `p`, `padj`,
#'   `method`.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 8,
#'                           group_close_ratio = list(all = 0.4),
#'                           n_types = 3, n_cells = 300, max_loc = 300,
#'                           test_type = "distribute", distance_param = 10,
#'                           seed = 2)
#' ps <- nhood_enrichment_per_sample(df, sample_key = "sample_id",
#'                                   group_key = "group",
#'                                   cluster_key = "cell_type",
#'                                   patient_key = "patient",
#'                                   neighbors.k = 10, n_perms = 30, n_jobs = 1)
#' crp <- setNames(rexp(8), unique(ps$patient))
#' head(associate_continuous(ps, crp, value = "log2_oe"))
associate_continuous <- function(per_sample_df, variable, value = "log2_oe",
                                 method = c("spearman", "lm", "lmm", "perm"),
                                 patient_key = "patient", covariates = NULL,
                                 pair_keys = NULL, n_perms = 2000, adjust = "BH",
                                 min_n = 4, seed = 1234) {
  method <- match.arg(method)
  d <- as.data.frame(per_sample_df)
  if (!value %in% colnames(d)) stop("value column '", value, "' not found.")
  if (!patient_key %in% colnames(d)) {
    if (!"sample_id" %in% colnames(d)) stop("need a '", patient_key, "' or 'sample_id' column.")
    d[[patient_key]] <- d$sample_id
  }
  if (is.character(variable) && length(variable) == 1L) {
    if (!variable %in% colnames(d)) stop("variable column '", variable, "' not found.")
    d$.x <- as.numeric(d[[variable]])
  } else {
    if (is.null(names(variable))) stop("`variable` must be a column name or a named numeric vector.")
    key <- if (all(unique(d[[patient_key]]) %in% names(variable))) d[[patient_key]] else d$sample_id
    d$.x <- as.numeric(variable[as.character(key)])
  }
  if (is.null(pair_keys)) pair_keys <- intersect(c("cluster_i", "cluster_j", "r"), colnames(d))
  if (method == "lmm" && !requireNamespace("lmerTest", quietly = TRUE)) {
    stop("method = 'lmm' requires the 'lmerTest' package.")
  }
  if (method == "perm") .local_seed(seed)
  d <- d[is.finite(d[[value]]) & is.finite(d$.x), , drop = FALSE]
  split_key <- interaction(d[, pair_keys, drop = FALSE], drop = TRUE, sep = "\r")
  one <- function(g) {
    pat <- g[[patient_key]]
    pm <- stats::aggregate(g[, c(value, ".x", covariates), drop = FALSE], by = list(patient = pat), FUN = mean)
    n <- nrow(pm)
    out <- data.frame(g[1, pair_keys, drop = FALSE], n_patients = n, estimate = NA_real_, p = NA_real_, row.names = NULL)
    if (n < min_n || stats::sd(pm$.x) == 0 || stats::sd(pm[[value]]) == 0) return(out)
    if (method %in% c("spearman", "perm")) {
      rho <- suppressWarnings(stats::cor(pm[[value]], pm$.x, method = "spearman"))
      out$estimate <- rho
      if (method == "spearman") {
        out$p <- suppressWarnings(stats::cor.test(pm[[value]], pm$.x, method = "spearman",
                                                  exact = n < 10 && !anyDuplicated(pm$.x) && !anyDuplicated(pm[[value]]))$p.value)
      } else {
        null <- replicate(n_perms, suppressWarnings(stats::cor(pm[[value]], sample(pm$.x), method = "spearman")))
        out$p <- (1 + sum(abs(null) >= abs(rho) - 1e-12)) / (n_perms + 1)
      }
    } else if (method == "lm" || max(table(pat)) < 2) {
      # lm on patient means (also used by "lmm" when every patient has one image)
      f <- stats::reformulate(c(".x", covariates), response = value)
      cf <- summary(stats::lm(f, data = pm))$coefficients
      if (".x" %in% rownames(cf)) { out$estimate <- cf[".x", 1]; out$p <- cf[".x", 4] }
    } else {
      f <- stats::as.formula(paste(value, "~", paste(c(".x", covariates), collapse = " + "), "+ (1 |", patient_key, ")"))
      fit <- tryCatch(suppressMessages(suppressWarnings(lmerTest::lmer(f, data = g))), error = function(e) NULL)
      if (!is.null(fit)) {
        cf <- summary(fit)$coefficients
        if (".x" %in% rownames(cf)) { out$estimate <- cf[".x", "Estimate"]; out$p <- cf[".x", "Pr(>|t|)"] }
      }
    }
    out
  }
  res <- do.call(rbind, lapply(split(d, split_key), one))
  res$padj <- stats::p.adjust(res$p, method = adjust)
  res$method <- method
  rownames(res) <- NULL
  res[order(res$p), , drop = FALSE]
}
