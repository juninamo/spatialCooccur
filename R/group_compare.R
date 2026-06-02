# ---- Disease / condition group comparison utilities ----
#
# These functions add a sample-aware layer on top of the single-object
# analyses in spatialCooccur_functions.R, so that scores can be computed
# per image (or aggregated per patient) and then compared between disease
# groups using standard or mixed-effects statistics, or a group-label
# permutation test.

#' @importFrom stats wilcox.test t.test p.adjust sd quantile pnorm aggregate as.formula lm rnorm
#' @importFrom utils head
NULL

utils::globalVariables(c(
  "cluster_i", "cluster_j", "value", "group", "patient", "sample_id",
  "cell_type"
))

# ---- Phase 0: input dispatch ---------------------------------------------

#' Build a sample design table for disease-group comparisons
#'
#' Construct a data.frame mapping image / sample identifiers to a disease
#' group label (and optionally a patient ID). This is the shared metadata
#' table consumed by the `*_per_sample` helpers and by `compare_groups()`.
#'
#' @param obj A Seurat object, a list of Seurat objects, or a data.frame.
#'   For Seurat input, `sample_key` is expected to be a `meta.data` column
#'   whose values match the image names in `obj@images`. For data.frame
#'   input, `sample_key` is a column of the data.frame.
#' @param sample_key Name of the column / image identifier.
#' @param group_key Name of the column carrying the disease group (or any
#'   condition label).
#' @param patient_key Optional column name for patient identifier, used as a
#'   random effect (or permutation block) downstream. NULL to omit.
#'
#' @return A data.frame with columns `sample_id`, `group`, `patient`, and
#'   `source_index` (which list element of the input the sample came from).
#' @export
build_sample_design <- function(obj, sample_key, group_key, patient_key = NULL) {
  if (inherits(obj, "Seurat")) {
    obj_list <- list(obj)
  } else if (is.list(obj) && length(obj) > 0L &&
             all(vapply(obj, inherits, logical(1), what = "Seurat"))) {
    obj_list <- obj
  } else if (is.data.frame(obj)) {
    df <- obj
    if (!sample_key %in% colnames(df)) {
      stop(sprintf("sample_key '%s' not found in data.frame columns", sample_key))
    }
    if (!group_key %in% colnames(df)) {
      stop(sprintf("group_key '%s' not found in data.frame columns", group_key))
    }
    sids <- unique(as.character(df[[sample_key]]))
    design <- data.frame(
      sample_id = sids,
      group = vapply(sids, function(s) {
        as.character(unique(df[[group_key]][df[[sample_key]] == s]))[1]
      }, character(1)),
      patient = if (is.null(patient_key)) NA_character_ else {
        if (!patient_key %in% colnames(df)) {
          stop(sprintf("patient_key '%s' not found in data.frame columns", patient_key))
        }
        vapply(sids, function(s) {
          as.character(unique(df[[patient_key]][df[[sample_key]] == s]))[1]
        }, character(1))
      },
      source_index = 1L,
      stringsAsFactors = FALSE
    )
    return(design)
  } else {
    stop("`obj` must be a Seurat object, a list of Seurat objects, or a data.frame.")
  }

  rows <- list()
  for (k in seq_along(obj_list)) {
    seu <- obj_list[[k]]
    md <- seu@meta.data
    if (!group_key %in% colnames(md)) {
      stop(sprintf("group_key '%s' not found in meta.data of object %d", group_key, k))
    }
    if (!is.null(patient_key) && !patient_key %in% colnames(md)) {
      stop(sprintf("patient_key '%s' not found in meta.data of object %d", patient_key, k))
    }
    for (img in names(seu@images)) {
      sel <- if (sample_key %in% colnames(md)) md[[sample_key]] == img else rep(TRUE, nrow(md))
      if (!any(sel)) next
      gv <- unique(as.character(md[[group_key]][sel]))
      if (length(gv) > 1L) {
        warning(sprintf(
          "Image '%s' (object %d) maps to multiple groups (%s); using first.",
          img, k, paste(gv, collapse = ", ")
        ))
      }
      pv <- if (is.null(patient_key)) NA_character_ else {
        v <- unique(as.character(md[[patient_key]][sel]))
        if (length(v) > 1L) {
          warning(sprintf(
            "Image '%s' (object %d) maps to multiple patients (%s); using first.",
            img, k, paste(v, collapse = ", ")
          ))
        }
        v[1]
      }
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = img, group = gv[1], patient = pv, source_index = k,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    stop("No samples were found. Check sample_key and the images slot of your object.")
  }
  do.call(rbind, rows)
}

# Internal: normalize any supported input into a per-sample list with
# `coords` (data.frame with x, y, cluster, rownames = cell IDs) plus
# `group` and `patient` metadata. Returns list(samples = ..., design = ...).
.as_sample_list <- function(obj, sample_key, cluster_key,
                            group_key = NULL, patient_key = NULL) {
  design <- if (!is.null(group_key)) {
    build_sample_design(obj, sample_key = sample_key,
                        group_key = group_key, patient_key = patient_key)
  } else {
    NULL
  }
  out <- list()

  if (is.data.frame(obj)) {
    df <- obj
    if (!cluster_key %in% colnames(df)) {
      stop(sprintf("cluster_key '%s' not found in data.frame columns", cluster_key))
    }
    if (!all(c("x", "y") %in% colnames(df))) {
      stop("data.frame input must contain 'x' and 'y' columns")
    }
    sids <- unique(as.character(df[[sample_key]]))
    for (s in sids) {
      sub <- df[as.character(df[[sample_key]]) == s, , drop = FALSE]
      cells <- if (!is.null(rownames(sub)) && !all(rownames(sub) == as.character(seq_len(nrow(sub))))) {
        rownames(sub)
      } else {
        paste0(s, "_cell", seq_len(nrow(sub)))
      }
      coords <- data.frame(
        x = sub$x, y = sub$y,
        cluster = as.character(sub[[cluster_key]]),
        row.names = cells, stringsAsFactors = FALSE
      )
      d <- if (!is.null(design)) design[design$sample_id == s, , drop = FALSE] else NULL
      out[[s]] <- list(
        coords = coords,
        group = if (!is.null(d) && nrow(d)) d$group[1] else NA_character_,
        patient = if (!is.null(d) && nrow(d)) d$patient[1] else NA_character_
      )
    }
    return(list(samples = out, design = design))
  }

  obj_list <- if (inherits(obj, "Seurat")) list(obj) else obj
  if (!is.list(obj_list)) {
    stop("`obj` must be a Seurat object, list of Seurat objects, or data.frame.")
  }

  for (k in seq_along(obj_list)) {
    seu <- obj_list[[k]]
    md <- seu@meta.data
    if (!cluster_key %in% colnames(md)) {
      stop(sprintf("cluster_key '%s' not found in meta.data of object %d", cluster_key, k))
    }
    for (img in names(seu@images)) {
      coords <- as.data.frame(seu[[img]]$centroids@coords)
      cells <- Seurat::Cells(seu[[img]])
      rownames(coords) <- cells
      sel <- if (sample_key %in% colnames(md)) md[[sample_key]] == img else rep(TRUE, nrow(md))
      md_sel <- md[sel, , drop = FALSE]
      common <- intersect(cells, rownames(md_sel))
      if (length(common) == 0L) {
        # fall back: assume order matches
        if (nrow(md_sel) != nrow(coords)) {
          stop(sprintf(
            "Image '%s' has %d cells but matching meta.data has %d rows.",
            img, nrow(coords), nrow(md_sel)
          ))
        }
        coords$cluster <- as.character(md_sel[[cluster_key]])
      } else {
        coords <- coords[common, , drop = FALSE]
        coords$cluster <- as.character(md_sel[common, cluster_key])
      }
      d <- if (!is.null(design)) {
        design[design$sample_id == img & design$source_index == k, , drop = FALSE]
      } else NULL
      out[[img]] <- list(
        coords = coords,
        group = if (!is.null(d) && nrow(d)) d$group[1] else NA_character_,
        patient = if (!is.null(d) && nrow(d)) d$patient[1] else NA_character_
      )
    }
  }
  list(samples = out, design = design)
}

# Internal: stack per-sample matrix results into a tidy long data.frame.
# `per_sample` is a named list of named lists; for each sample we take
# the matrices specified by `value_map`, vectorize them, and stack.
# value_map: named character vector list-element-name -> output-column-name.
.stack_pair_results <- function(per_sample, design, value_map) {
  rows <- list()
  for (sid in names(per_sample)) {
    res <- per_sample[[sid]]
    if (is.null(res)) next
    mat0 <- res[[names(value_map)[1]]]
    if (is.null(mat0)) next
    rn <- rownames(mat0); cn <- colnames(mat0)
    grid <- expand.grid(cluster_i = rn, cluster_j = cn, stringsAsFactors = FALSE)
    for (k in seq_along(value_map)) {
      src <- names(value_map)[k]
      dst <- unname(value_map[k])
      grid[[dst]] <- as.vector(res[[src]])
    }
    grid$sample_id <- sid
    rows[[sid]] <- grid
  }
  if (length(rows) == 0L) return(NULL)
  out <- do.call(rbind, rows)
  out$cluster_i <- sub("^Cluster", "", out$cluster_i)
  out$cluster_j <- sub("^Cluster", "", out$cluster_j)
  if (!is.null(design)) {
    out <- merge(out, design[, c("sample_id", "group", "patient")],
                 by = "sample_id", all.x = TRUE)
  }
  rownames(out) <- NULL
  out
}

# Internal: aggregate sample-level rows to patient-level rows by mean.
.aggregate_to_patient <- function(tidy_df) {
  if (!"patient" %in% colnames(tidy_df) || all(is.na(tidy_df$patient))) {
    stop("patient column is missing or all NA; cannot aggregate to patient level.")
  }
  num_cols <- vapply(tidy_df, is.numeric, logical(1))
  num_cols <- names(num_cols)[num_cols]
  by_cols <- list(
    patient = tidy_df$patient,
    cluster_i = if ("cluster_i" %in% colnames(tidy_df)) tidy_df$cluster_i else rep("", nrow(tidy_df)),
    cluster_j = if ("cluster_j" %in% colnames(tidy_df)) tidy_df$cluster_j else rep("", nrow(tidy_df)),
    group = tidy_df$group
  )
  agg <- aggregate(tidy_df[, num_cols, drop = FALSE], by = by_cols,
                   FUN = function(x) mean(x, na.rm = TRUE))
  agg$sample_id <- agg$patient
  agg
}

# ---- Phase 1a: nhood_enrichment_per_sample -------------------------------

#' Per-sample neighborhood enrichment
#'
#' Run [nhood_enrichment()] independently for each sample (image) in the
#' input, returning a tidy long data.frame with one row per
#' `sample_id x cluster_i x cluster_j`. The cluster factor levels are
#' harmonized across samples so the resulting table is suitable for
#' downstream group comparison with [compare_groups()].
#'
#' @inheritParams nhood_enrichment
#' @param obj A Seurat object, a list of Seurat objects, or a data.frame
#'   with x, y, cluster_key, and sample_key columns.
#' @param sample_key Column / image identifier defining a sample.
#' @param group_key Column carrying the disease / condition label.
#' @param patient_key Optional column with patient ID (random effect /
#'   permutation block).
#' @param unit "image" returns one row per image; "patient" averages across
#'   images within a patient (requires `patient_key`).
#' @param cluster_levels Optional character vector of cluster levels to use
#'   as common dimnames across samples. Defaults to the union across all
#'   samples.
#'
#' @return A data.frame (also tagged with class `spatialCooccurSample`) with
#'   columns `sample_id`, `cluster_i`, `cluster_j`, `zscore`, `count`,
#'   `group`, `patient`.
#' @export
nhood_enrichment_per_sample <- function(obj, sample_key, group_key, cluster_key,
                                        patient_key = NULL,
                                        unit = c("image", "patient"),
                                        cluster_levels = NULL,
                                        neighbors.k = 30,
                                        connectivity_key = "nn",
                                        transformation = TRUE,
                                        n_perms = 100,
                                        seed = 1938493,
                                        n_jobs = 1) {
  unit <- match.arg(unit)
  sl <- .as_sample_list(obj, sample_key = sample_key,
                        cluster_key = cluster_key,
                        group_key = group_key, patient_key = patient_key)
  samples <- sl$samples
  design <- sl$design

  if (is.null(cluster_levels)) {
    cluster_levels <- sort(unique(unlist(lapply(samples, function(s) unique(s$coords$cluster)))))
  }

  per_sample <- list()
  for (sid in names(samples)) {
    coords <- samples[[sid]]$coords
    if (nrow(coords) < 2L) {
      warning(sprintf("Sample '%s' has fewer than 2 cells; skipping.", sid))
      next
    }
    df_in <- data.frame(
      x = coords$x, y = coords$y,
      stringsAsFactors = FALSE,
      row.names = rownames(coords)
    )
    df_in[[cluster_key]] <- factor(coords$cluster, levels = cluster_levels)
    res <- tryCatch(
      nhood_enrichment(df_in, cluster_key = cluster_key,
                       neighbors.k = neighbors.k,
                       connectivity_key = connectivity_key,
                       transformation = transformation,
                       n_perms = n_perms,
                       seed = seed,
                       n_jobs = n_jobs),
      error = function(e) {
        warning(sprintf("nhood_enrichment failed for sample '%s': %s", sid, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) per_sample[[sid]] <- res
  }

  tidy_df <- .stack_pair_results(per_sample, design,
                                 value_map = c(zscore = "zscore", count = "count"))
  if (is.null(tidy_df) || nrow(tidy_df) == 0L) {
    stop("No per-sample results were produced.")
  }
  if (unit == "patient") tidy_df <- .aggregate_to_patient(tidy_df)
  attr(tidy_df, "spatial_design") <- design
  attr(tidy_df, "value_columns") <- c("zscore", "count")
  class(tidy_df) <- c("spatialCooccurSample", class(tidy_df))
  tidy_df
}

# ---- Phase 1b: cooccur_local_per_sample ---------------------------------

#' Per-sample local co-occurrence summary
#'
#' Run [cooccur_local()] for one (cluster_x, cluster_y) pair on each sample
#' and summarize the per-cell scores at the sample level. Output is a tidy
#' data.frame suitable for group comparison.
#'
#' @inheritParams cooccur_local
#' @param obj A Seurat object, a list of Seurat objects, or a data.frame.
#' @param sample_key,group_key,cluster_key,patient_key Sample / group /
#'   cluster / patient column names. See [build_sample_design()].
#' @param unit "image" or "patient" (averaged across images of the same
#'   patient).
#' @param summarize Character vector of summary statistics to compute:
#'   "mean", "q90" (90th percentile), and / or "pos_rate" (fraction of
#'   cells with score > 0).
#'
#' @return A data.frame with one row per sample carrying the requested
#'   summary statistics.
#' @export
cooccur_local_per_sample <- function(obj, sample_key, group_key, cluster_key,
                                     cluster_x, cluster_y,
                                     patient_key = NULL,
                                     unit = c("image", "patient"),
                                     neighbors.k = 20, radius = 30, maxnsteps = 1,
                                     summarize = c("mean", "q90", "pos_rate")) {
  unit <- match.arg(unit)
  summarize <- match.arg(summarize, several.ok = TRUE)
  sl <- .as_sample_list(obj, sample_key = sample_key,
                        cluster_key = cluster_key,
                        group_key = group_key, patient_key = patient_key)
  samples <- sl$samples
  design <- sl$design

  rows <- list()
  for (sid in names(samples)) {
    coords <- samples[[sid]]$coords
    if (nrow(coords) < 2L) next
    df_in <- data.frame(
      x = coords$x, y = coords$y,
      cell_type = as.character(coords$cluster),
      row.names = rownames(coords),
      stringsAsFactors = FALSE
    )
    res <- tryCatch(
      cooccur_local(df_in, cluster_x = cluster_x, cluster_y = cluster_y,
                    neighbors.k = neighbors.k, radius = radius,
                    maxnsteps = maxnsteps),
      error = function(e) {
        warning(sprintf("cooccur_local failed for sample '%s': %s", sid, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) next
    sc <- as.numeric(res[[1]])
    row <- data.frame(
      sample_id = sid,
      cluster_i = cluster_x,
      cluster_j = cluster_y,
      stringsAsFactors = FALSE
    )
    if ("mean" %in% summarize) row$mean <- mean(sc, na.rm = TRUE)
    if ("q90" %in% summarize) row$q90 <- as.numeric(quantile(sc, 0.9, na.rm = TRUE))
    if ("pos_rate" %in% summarize) row$pos_rate <- mean(sc > 0, na.rm = TRUE)
    rows[[sid]] <- row
  }
  if (length(rows) == 0L) stop("No per-sample results were produced.")
  tidy_df <- do.call(rbind, rows)
  if (!is.null(design)) {
    tidy_df <- merge(tidy_df, design[, c("sample_id", "group", "patient")],
                     by = "sample_id", all.x = TRUE)
  }
  if (unit == "patient") tidy_df <- .aggregate_to_patient(tidy_df)
  attr(tidy_df, "spatial_design") <- design
  attr(tidy_df, "value_columns") <- intersect(summarize, colnames(tidy_df))
  class(tidy_df) <- c("spatialCooccurSample", class(tidy_df))
  tidy_df
}

# ---- Phase 1b: cooccur_ratio_per_sample ---------------------------------

# Internal: radius-based co-occurrence count + ratio for a single sample.
.radius_count_one <- function(coords_xy, clusters, all_clusters, radius, k = 30) {
  res <- RANN::nn2(data = coords_xy, query = coords_xy,
                   searchtype = "radius", radius = radius, k = k)
  co_occur_count <- matrix(
    0,
    nrow = length(all_clusters), ncol = length(all_clusters),
    dimnames = list(paste0("Cluster", all_clusters),
                    paste0("Cluster", all_clusters))
  )
  for (i in seq_len(nrow(coords_xy))) {
    nb <- res$nn.idx[i, ]
    nb <- nb[nb != i & nb > 0]
    if (!length(nb)) next
    ci <- paste0("Cluster", clusters[i])
    tab <- table(paste0("Cluster", clusters[nb]))
    for (cn in names(tab)) {
      co_occur_count[ci, cn] <- co_occur_count[ci, cn] + tab[[cn]]
    }
  }
  ratio_mat <- compute_co_occurrence_ratio(co_occur_count)
  list(co_occur_count = co_occur_count, ratio_mat = ratio_mat)
}

#' Per-sample radius-based co-occurrence ratio
#'
#' Compute the radius-based co-occurrence count and enrichment ratio per
#' sample, with cluster dimnames harmonized across samples.
#'
#' @inheritParams calc_co_occurrence_for_radius
#' @param obj A Seurat object, list of Seurat objects, or data.frame.
#' @param sample_key,group_key,cluster_key,patient_key See
#'   [build_sample_design()].
#' @param unit "image" or "patient".
#' @param cluster_levels Optional vector of cluster levels.
#'
#' @return A data.frame with one row per `sample_id x cluster_i x cluster_j`,
#'   columns `ratio`, `count`, `group`, `patient`.
#' @export
cooccur_ratio_per_sample <- function(obj, sample_key, group_key, cluster_key,
                                     patient_key = NULL,
                                     unit = c("image", "patient"),
                                     radius = 30, k = 30,
                                     cluster_levels = NULL) {
  unit <- match.arg(unit)
  sl <- .as_sample_list(obj, sample_key = sample_key,
                        cluster_key = cluster_key,
                        group_key = group_key, patient_key = patient_key)
  samples <- sl$samples
  design <- sl$design

  if (is.null(cluster_levels)) {
    cluster_levels <- sort(unique(unlist(lapply(samples, function(s) unique(s$coords$cluster)))))
  }

  per_sample <- list()
  for (sid in names(samples)) {
    coords <- samples[[sid]]$coords
    if (nrow(coords) < 2L) next
    res <- tryCatch(
      .radius_count_one(coords_xy = as.matrix(coords[, c("x", "y")]),
                        clusters = as.character(coords$cluster),
                        all_clusters = cluster_levels,
                        radius = radius, k = k),
      error = function(e) {
        warning(sprintf("ratio computation failed for sample '%s': %s",
                        sid, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) per_sample[[sid]] <- res
  }

  tidy_df <- .stack_pair_results(per_sample, design,
                                 value_map = c(ratio_mat = "ratio",
                                               co_occur_count = "count"))
  if (is.null(tidy_df) || nrow(tidy_df) == 0L) {
    stop("No per-sample results were produced.")
  }
  if (unit == "patient") tidy_df <- .aggregate_to_patient(tidy_df)
  attr(tidy_df, "spatial_design") <- design
  attr(tidy_df, "value_columns") <- c("ratio", "count")
  class(tidy_df) <- c("spatialCooccurSample", class(tidy_df))
  tidy_df
}

# ---- Phase 1b: interaction_spot_per_sample ------------------------------

#' Per-sample interaction-spot summary
#'
#' Apply [search_interaction_spot()] to each image of a Seurat object and
#' summarize how many connected interaction spots were detected and their
#' mean size, ready for group comparison.
#'
#' @param seurat_object A Seurat object (lists are not supported here).
#' @param sample_key,group_key,patient_key Column names; see
#'   [build_sample_design()].
#' @param cluster_col Cluster column in meta.data.
#' @param target_cluster Target cluster(s) of interest.
#' @param cell_id Vector of cell IDs to include. Defaults to all cells.
#' @param radius Radius defining neighborhood.
#' @param n_min Minimum number of cells per spot.
#' @param neighbors.k Max neighbors to consider.
#'
#' @return A data.frame with one row per sample carrying `n_spots` and
#'   `mean_spot_size` (number of cells).
#' @export
interaction_spot_per_sample <- function(seurat_object, sample_key, group_key,
                                        cluster_col, target_cluster,
                                        cell_id = NULL,
                                        patient_key = NULL,
                                        radius, n_min, neighbors.k = 200) {
  if (!inherits(seurat_object, "Seurat")) {
    stop("interaction_spot_per_sample currently supports Seurat object input only.")
  }
  design <- build_sample_design(seurat_object, sample_key = sample_key,
                                group_key = group_key, patient_key = patient_key)
  if (is.null(cell_id)) {
    if ("cell" %in% colnames(seurat_object@meta.data)) {
      cell_id <- as.character(seurat_object@meta.data$cell)
    } else {
      cell_id <- rownames(seurat_object@meta.data)
    }
  }
  rows <- list()
  for (fov in names(seurat_object@images)) {
    spots <- tryCatch(
      search_interaction_spot(seurat_object, fov = fov, radius = radius, n_min = n_min,
                              neighbors.k = neighbors.k, cell_id = cell_id,
                              cluster_col = cluster_col,
                              target_cluster = target_cluster),
      error = function(e) {
        warning(sprintf("search_interaction_spot failed for '%s': %s",
                        fov, conditionMessage(e)))
        NULL
      }
    )
    n_spots <- if (!is.null(spots)) length(unique(spots$cluster_id)) else 0L
    mean_size <- if (!is.null(spots) && nrow(spots) > 0L) {
      sizes <- unique(spots[, c("cluster_id", "n_all_cells")])
      mean(sizes$n_all_cells, na.rm = TRUE)
    } else NA_real_
    d <- design[design$sample_id == fov, , drop = FALSE]
    rows[[fov]] <- data.frame(
      sample_id = fov,
      group = if (nrow(d)) d$group[1] else NA_character_,
      patient = if (nrow(d)) d$patient[1] else NA_character_,
      target_cluster = paste(target_cluster, collapse = ","),
      n_spots = n_spots,
      mean_spot_size = mean_size,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  attr(out, "spatial_design") <- design
  attr(out, "value_columns") <- c("n_spots", "mean_spot_size")
  class(out) <- c("spatialCooccurSample", class(out))
  out
}

# ---- Phase 2: compare_groups --------------------------------------------

#' Compare disease groups across samples
#'
#' Given a tidy per-sample data.frame produced by one of the
#' `*_per_sample()` helpers, run a per-cluster-pair statistical test
#' between groups and return tidy results with multiple-testing adjustment.
#'
#' @param per_sample_df Tidy data.frame, typically the output of
#'   [nhood_enrichment_per_sample()], [cooccur_ratio_per_sample()],
#'   [cooccur_local_per_sample()], or [interaction_spot_per_sample()].
#' @param value Name of the column to test (e.g. "zscore", "ratio", "mean",
#'   "n_spots").
#' @param group_key Name of the group column. Defaults to "group".
#' @param patient_key Optional patient column. Required for `method = "lmm"`
#'   to use as a random effect; used as the permutation block for
#'   `method = "perm"`.
#' @param method Statistical test:
#'   * "wilcox" — Wilcoxon rank-sum (two groups)
#'   * "t" — Welch's t-test (two groups)
#'   * "lmm" — linear mixed model `value ~ group + (1|patient)` via
#'     `lme4::lmer`, with Wald-z p-values. Requires the `lme4` package.
#'   * "perm" — group-label permutation test on the mean difference
#'     (blocked by patient if `patient_key` is supplied).
#' @param n_perms Number of permutations for `method = "perm"`.
#' @param adjust Multiple-testing adjustment method passed to
#'   [stats::p.adjust()].
#' @param pair_keys Column names that together identify a cluster pair.
#'   Defaults to `c("cluster_i", "cluster_j")`. Set to a single column name
#'   for cases like interaction_spot_per_sample (e.g. "target_cluster").
#' @param ref_group Optional name of the reference group. If supplied, the
#'   effect is reported as `mean_other - mean_ref` (positive = higher in
#'   the non-reference group). When `NULL` (default), groups are sorted
#'   alphabetically and the second is treated as the test group.
#' @param seed Random seed for the permutation test.
#'
#' @return A data.frame with the cluster pair columns, group means, effect
#'   size (test group minus reference group), test statistic, raw p-value,
#'   and adjusted p-value.
#' @export
compare_groups <- function(per_sample_df,
                           value = "zscore",
                           group_key = "group",
                           patient_key = NULL,
                           method = c("wilcox", "t", "lmm", "perm"),
                           n_perms = 1000,
                           adjust = "BH",
                           pair_keys = c("cluster_i", "cluster_j"),
                           ref_group = NULL,
                           seed = 1234) {
  method <- match.arg(method)
  if (!value %in% colnames(per_sample_df)) {
    stop(sprintf("value column '%s' not in data.frame", value))
  }
  if (!group_key %in% colnames(per_sample_df)) {
    stop(sprintf("group_key column '%s' not in data.frame", group_key))
  }
  if (!all(pair_keys %in% colnames(per_sample_df))) {
    missing_keys <- setdiff(pair_keys, colnames(per_sample_df))
    stop(sprintf("pair_keys not all present: %s", paste(missing_keys, collapse = ", ")))
  }
  if (method == "lmm" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("method = 'lmm' requires the 'lme4' package; install it or choose another method.")
  }
  if (method == "lmm" && is.null(patient_key)) {
    warning("method = 'lmm' without patient_key is equivalent to OLS; passing patient_key is recommended.")
  }

  set.seed(seed)
  df <- as.data.frame(per_sample_df)

  groups_present <- unique(as.character(df[[group_key]]))
  groups_present <- groups_present[!is.na(groups_present)]
  if (length(groups_present) < 2L) {
    stop("Need at least 2 groups to compare. Found: ", paste(groups_present, collapse = ", "))
  }
  # Determine reference (g1) and test (g2) levels.
  if (!is.null(ref_group)) {
    if (!ref_group %in% groups_present) {
      stop(sprintf("ref_group '%s' not found in '%s'. Available: %s",
                   ref_group, group_key, paste(groups_present, collapse = ", ")))
    }
    g1 <- ref_group
    other <- setdiff(groups_present, ref_group)
    if (length(other) > 1L && method %in% c("wilcox", "t", "perm")) {
      warning(sprintf(
        "Method '%s' is for 2-group comparison; using ref '%s' vs '%s'.",
        method, g1, other[1]
      ))
    }
    g2 <- other[1]
  } else {
    groups_sorted <- sort(groups_present)
    if (length(groups_sorted) > 2L && method %in% c("wilcox", "t", "perm")) {
      warning(sprintf(
        "Method '%s' is for 2-group comparison; using first 2 sorted levels (%s vs %s). Pass ref_group to control this.",
        method, groups_sorted[1], groups_sorted[2]
      ))
    }
    g1 <- groups_sorted[1]
    g2 <- groups_sorted[2]
  }

  # Unique pairs
  pair_ids <- unique(df[, pair_keys, drop = FALSE])
  rownames(pair_ids) <- NULL

  results <- vector("list", nrow(pair_ids))
  for (i in seq_len(nrow(pair_ids))) {
    pi <- pair_ids[i, , drop = FALSE]
    sub <- merge(df, pi, by = pair_keys)
    sub <- sub[is.finite(sub[[value]]), , drop = FALSE]
    sub <- sub[as.character(sub[[group_key]]) %in% c(g1, g2), , drop = FALSE]
    if (nrow(sub) < 2L) next
    g <- as.character(sub[[group_key]])
    if (length(unique(g)) < 2L) next

    res_row <- as.list(pi)
    res_row$n_total <- nrow(sub)
    res_row[[paste0("n_", g1)]] <- sum(g == g1)
    res_row[[paste0("n_", g2)]] <- sum(g == g2)
    v <- sub[[value]]
    m1 <- mean(v[g == g1], na.rm = TRUE)
    m2 <- mean(v[g == g2], na.rm = TRUE)
    res_row[[paste0("mean_", g1)]] <- m1
    res_row[[paste0("mean_", g2)]] <- m2
    res_row$effect <- m2 - m1

    test <- list(stat = NA_real_, p = NA_real_)
    if (method == "wilcox") {
      tt <- tryCatch(
        wilcox.test(v ~ factor(g, levels = c(g1, g2)), exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(tt)) { test$stat <- unname(tt$statistic); test$p <- tt$p.value }
    } else if (method == "t") {
      tt <- tryCatch(t.test(v ~ factor(g, levels = c(g1, g2))), error = function(e) NULL)
      if (!is.null(tt)) { test$stat <- unname(tt$statistic); test$p <- tt$p.value }
    } else if (method == "lmm") {
      sub$.g <- factor(g, levels = c(g1, g2))
      use_lmm <- !is.null(patient_key) && patient_key %in% colnames(sub) &&
        length(unique(sub[[patient_key]])) < nrow(sub)
      if (use_lmm) {
        sub$.p <- sub[[patient_key]]
        fit <- tryCatch(
          lme4::lmer(as.formula(sprintf("%s ~ .g + (1 | .p)", value)),
                     data = sub, REML = TRUE),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          cf <- summary(fit)$coefficients
          if (nrow(cf) >= 2L) {
            zval <- cf[2, "t value"]
            test$stat <- zval
            test$p <- 2 * pnorm(-abs(zval))
          }
        }
      } else {
        fit <- tryCatch(lm(as.formula(sprintf("%s ~ .g", value)), data = sub),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          cf <- summary(fit)$coefficients
          if (nrow(cf) >= 2L) {
            test$stat <- cf[2, "t value"]
            test$p <- cf[2, "Pr(>|t|)"]
          }
        }
      }
    } else if (method == "perm") {
      observed <- m2 - m1
      if (!is.null(patient_key) && patient_key %in% colnames(sub)) {
        pat <- as.character(sub[[patient_key]])
        # Map each patient to its group, shuffle group labels at patient level
        pat_levels <- unique(pat)
        pat_to_group <- vapply(pat_levels, function(p) g[which(pat == p)[1]], character(1))
        names(pat_to_group) <- pat_levels
        null_diffs <- replicate(n_perms, {
          shuffled <- sample(pat_to_group)
          names(shuffled) <- pat_levels
          new_g <- shuffled[pat]
          mean(v[new_g == g2], na.rm = TRUE) - mean(v[new_g == g1], na.rm = TRUE)
        })
      } else {
        null_diffs <- replicate(n_perms, {
          new_g <- sample(g)
          mean(v[new_g == g2], na.rm = TRUE) - mean(v[new_g == g1], na.rm = TRUE)
        })
      }
      test$stat <- observed
      test$p <- (sum(abs(null_diffs) >= abs(observed), na.rm = TRUE) + 1) /
        (sum(!is.na(null_diffs)) + 1)
    }

    res_row$statistic <- test$stat
    res_row$p <- test$p
    results[[i]] <- as.data.frame(res_row, stringsAsFactors = FALSE)
  }
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0L) {
    warning("No comparable cluster pairs found.")
    return(invisible(NULL))
  }
  out <- do.call(rbind, results)
  out$padj <- p.adjust(out$p, method = adjust)
  out <- out[order(out$p), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "method") <- method
  attr(out, "groups") <- c(g1, g2)
  attr(out, "value") <- value
  out
}

# ---- Phase 3: multi-group simulator -------------------------------------

#' Generate multi-sample simulated data with disease-group structure
#'
#' Wraps [generate_sim()] to produce N samples per group, each with a
#' (possibly noised) group-specific `close_ratio`. Returns a single tidy
#' data.frame with `x`, `y`, `cell_type`, `sample_id`, `group`, `patient`
#' columns — directly consumable by [nhood_enrichment_per_sample()] and the
#' other `*_per_sample` helpers.
#'
#' @param n_samples_per_group Integer, number of samples to generate per
#'   group.
#' @param group_close_ratio Named list of base `close_ratio` values, one
#'   entry per group, e.g. `list(disease = 0.8, control = 0.2)`.
#' @param n_types,max_loc,n_cells,test_type,distance_param Passed through
#'   to [generate_sim()].
#' @param between_sample_noise SD of Gaussian noise added to the
#'   per-sample `close_ratio` around the group baseline (clipped to [0,1]).
#' @param seed Random seed (controls both the noise and the per-sample
#'   seeds passed to `generate_sim`).
#'
#' @return A data.frame.
#' @export
generate_sim_groups <- function(n_samples_per_group = 3,
                                group_close_ratio = list(disease = 0.8, control = 0.2),
                                n_types = 10,
                                max_loc = 800,
                                n_cells = 1500,
                                test_type = "circle",
                                distance_param = 50,
                                between_sample_noise = 0.05,
                                seed = 1234) {
  set.seed(seed)
  groups <- names(group_close_ratio)
  if (is.null(groups) || any(groups == "")) {
    stop("group_close_ratio must be a named list, e.g. list(disease = 0.8, control = 0.2)")
  }
  out <- list()
  for (g in groups) {
    base_ratio <- group_close_ratio[[g]]
    for (s in seq_len(n_samples_per_group)) {
      sid <- paste0(g, "_", s)
      ratio_s <- base_ratio + rnorm(1, mean = 0, sd = between_sample_noise)
      ratio_s <- max(0, min(1, ratio_s))
      df <- generate_sim(close_ratio = ratio_s,
                         n_types = n_types,
                         max_loc = max_loc,
                         n_cells = n_cells,
                         test_type = test_type,
                         distance_param = distance_param,
                         seed = sample.int(.Machine$integer.max, 1))
      df$sample_id <- sid
      df$group <- g
      df$patient <- sid
      rownames(df) <- paste0(sid, "_cell", seq_len(nrow(df)))
      out[[sid]] <- df
    }
  }
  df_all <- do.call(rbind, out)
  rownames(df_all) <- unlist(lapply(out, rownames))
  df_all
}
