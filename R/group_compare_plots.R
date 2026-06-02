# ---- Disease / condition group comparison plots ----
#
# Lightweight ggplot2 helpers for visualizing the output of compare_groups()
# and the per-sample tidy data.frames. ggplot2 is a Suggests dependency:
# each function checks for availability and errors out with a helpful
# message if it is missing.

utils::globalVariables(c(
  "effect", "padj", "neg_log10_p", "sig", "score", ".pair"
))

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This function requires the 'ggplot2' package. Install it with install.packages('ggplot2').")
  }
}

#' Heatmap of the disease-group effect for every cluster pair
#'
#' Plot `effect` (or another column from a [compare_groups()] result) as a
#' cluster_i x cluster_j heatmap, with optional significance markers.
#'
#' @param compare_df Output of [compare_groups()].
#' @param value Column to map to fill. Defaults to "effect".
#' @param significance Column name carrying p-values to overlay as labels
#'   (set to NULL to skip).
#' @param sig_threshold Adjusted p-value threshold for the asterisk marker.
#' @param pair_keys Two column names identifying the cluster pair on the
#'   x and y axes. Defaults to `c("cluster_i", "cluster_j")`.
#' @param palette One of "RdBu" (diverging, default) or "viridis".
#' @param limits Optional length-2 numeric vector for the fill scale limits.
#'   If NULL, symmetric limits around 0 are used for diverging palettes.
#'
#' @return A ggplot object.
#' @export
plot_group_delta_heatmap <- function(compare_df,
                                     value = "effect",
                                     significance = "padj",
                                     sig_threshold = 0.05,
                                     pair_keys = c("cluster_i", "cluster_j"),
                                     palette = c("RdBu", "viridis"),
                                     limits = NULL) {
  .require_ggplot2()
  palette <- match.arg(palette)
  if (!all(pair_keys %in% colnames(compare_df))) {
    stop("pair_keys not all present in compare_df.")
  }
  if (!value %in% colnames(compare_df)) {
    stop(sprintf("value '%s' not in compare_df.", value))
  }

  df <- as.data.frame(compare_df)
  df$.x <- df[[pair_keys[1]]]
  df$.y <- df[[pair_keys[2]]]
  df$.v <- df[[value]]

  groups <- attr(compare_df, "groups")
  subtitle <- if (!is.null(groups)) {
    sprintf("%s: %s vs %s (ref)", value, groups[2], groups[1])
  } else {
    value
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.x, y = .data$.y, fill = .data$.v)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = pair_keys[1], y = pair_keys[2], fill = value,
                  title = "Disease-group effect heatmap",
                  subtitle = subtitle) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (palette == "RdBu") {
    lim <- if (!is.null(limits)) limits else {
      mx <- max(abs(df$.v), na.rm = TRUE)
      if (!is.finite(mx) || mx == 0) c(-1, 1) else c(-mx, mx)
    }
    p <- p + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                           high = "#B2182B", midpoint = 0,
                                           limits = lim)
  } else {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      p <- p + ggplot2::scale_fill_viridis_c(limits = limits)
    } else {
      p <- p + ggplot2::scale_fill_gradient(low = "#FDE725", high = "#440154",
                                            limits = limits)
    }
  }

  if (!is.null(significance) && significance %in% colnames(df)) {
    df$.sig <- ifelse(!is.na(df[[significance]]) & df[[significance]] < sig_threshold, "*", "")
    p <- p + ggplot2::geom_text(data = df, ggplot2::aes(label = .data$.sig),
                                size = 5, vjust = 0.75)
  }
  p
}

#' Per-sample boxplot of one (or several) cluster pair(s) across groups
#'
#' @param per_sample_df Tidy per-sample data.frame from a `*_per_sample()`
#'   helper.
#' @param value Numeric column to plot. Defaults to "zscore".
#' @param group_key Column with disease group. Defaults to "group".
#' @param pair_keys Column names identifying the cluster pair.
#' @param pairs Optional data.frame (or list of length-2 vectors) selecting
#'   which cluster pairs to plot. If NULL, plots all pairs as a facetted
#'   grid (use sparingly for many pairs).
#' @param point Logical, overlay individual sample points.
#' @param add_p Logical, overlay a Wilcoxon p-value on each panel.
#' @param ref_group Optional ref group for the displayed p-value sign /
#'   ordering. Passed to [compare_groups()] when computing p-values.
#'
#' @return A ggplot object.
#' @export
plot_pair_boxplot <- function(per_sample_df,
                              value = "zscore",
                              group_key = "group",
                              pair_keys = c("cluster_i", "cluster_j"),
                              pairs = NULL,
                              point = TRUE,
                              add_p = FALSE,
                              ref_group = NULL) {
  .require_ggplot2()
  df <- as.data.frame(per_sample_df)
  if (!value %in% colnames(df)) stop(sprintf("value '%s' not in data.frame.", value))
  if (!group_key %in% colnames(df)) stop(sprintf("group_key '%s' not in data.frame.", group_key))
  if (!all(pair_keys %in% colnames(df))) stop("pair_keys not all present.")

  if (!is.null(pairs)) {
    if (is.list(pairs) && !is.data.frame(pairs)) {
      pairs <- do.call(rbind, lapply(pairs, function(p) {
        setNames(as.data.frame(t(p), stringsAsFactors = FALSE), pair_keys)
      }))
    }
    df <- merge(df, pairs, by = pair_keys)
    if (nrow(df) == 0L) stop("No rows match the requested cluster pairs.")
  }

  df$.pair <- if (length(pair_keys) == 1L) {
    as.character(df[[pair_keys]])
  } else {
    paste(df[[pair_keys[1]]], df[[pair_keys[2]]], sep = " - ")
  }
  df$.score <- df[[value]]
  df$.group <- factor(df[[group_key]])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.group, y = .data$.score, fill = .data$.group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::labs(x = group_key, y = value, fill = group_key) +
    ggplot2::theme_minimal()
  if (point) {
    p <- p + ggplot2::geom_jitter(width = 0.15, alpha = 0.7, size = 1.4)
  }
  if (length(unique(df$.pair)) > 1L) {
    p <- p + ggplot2::facet_wrap(~ .pair, scales = "free_y")
  } else {
    p <- p + ggplot2::ggtitle(unique(df$.pair))
  }

  if (add_p) {
    cmp <- tryCatch(
      compare_groups(df, value = value, group_key = group_key,
                     method = "wilcox", adjust = "none",
                     pair_keys = pair_keys, ref_group = ref_group),
      error = function(e) NULL
    )
    if (!is.null(cmp) && nrow(cmp) > 0L) {
      cmp$.pair <- if (length(pair_keys) == 1L) {
        as.character(cmp[[pair_keys]])
      } else {
        paste(cmp[[pair_keys[1]]], cmp[[pair_keys[2]]], sep = " - ")
      }
      cmp$.lab <- sprintf("p = %.2g", cmp$p)
      anno <- merge(
        aggregate(.score ~ .pair, data = df, FUN = function(x) max(x, na.rm = TRUE)),
        cmp[, c(".pair", ".lab")], by = ".pair"
      )
      anno$.group <- factor(levels(df$.group)[1], levels = levels(df$.group))
      p <- p + ggplot2::geom_text(data = anno,
                                  ggplot2::aes(x = .data$.group, y = .data$.score,
                                               label = .data$.lab),
                                  hjust = 0, vjust = -0.2, size = 3, inherit.aes = FALSE)
    }
  }
  p
}

#' Volcano plot of a compare_groups result
#'
#' @param compare_df Output of [compare_groups()].
#' @param effect_col Column for the x axis. Defaults to "effect".
#' @param p_col Column for the y axis (will be -log10 transformed).
#'   Defaults to "padj".
#' @param sig_threshold Adjusted p-value threshold for highlighting.
#' @param label_top Integer, number of top pairs to label by p-value.
#'   Set to 0 to skip labels.
#' @param pair_keys Columns identifying the cluster pair (used for labels).
#'
#' @return A ggplot object.
#' @export
plot_volcano_groups <- function(compare_df,
                                effect_col = "effect",
                                p_col = "padj",
                                sig_threshold = 0.05,
                                label_top = 10,
                                pair_keys = c("cluster_i", "cluster_j")) {
  .require_ggplot2()
  if (!effect_col %in% colnames(compare_df)) {
    stop(sprintf("effect_col '%s' not in compare_df.", effect_col))
  }
  if (!p_col %in% colnames(compare_df)) {
    stop(sprintf("p_col '%s' not in compare_df.", p_col))
  }
  df <- as.data.frame(compare_df)
  df$.effect <- df[[effect_col]]
  df$.p <- df[[p_col]]
  df$.neg_log10 <- -log10(pmax(df$.p, .Machine$double.eps))
  df$.sig <- factor(
    ifelse(is.na(df$.p), "NA",
           ifelse(df$.p < sig_threshold,
                  ifelse(df$.effect > 0, "up", "down"),
                  "ns")),
    levels = c("up", "down", "ns", "NA")
  )

  cols <- c(up = "#B2182B", down = "#2166AC", ns = "grey70", `NA` = "grey90")

  groups <- attr(compare_df, "groups")
  subtitle <- if (!is.null(groups)) {
    sprintf("effect = mean(%s) - mean(%s)", groups[2], groups[1])
  } else {
    NULL
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.effect, y = .data$.neg_log10, color = .data$.sig)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.6) +
    ggplot2::geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = effect_col, y = sprintf("-log10(%s)", p_col),
                  color = sprintf("%s < %.2g", p_col, sig_threshold),
                  title = "Disease-group volcano plot",
                  subtitle = subtitle) +
    ggplot2::theme_minimal()

  if (label_top > 0L && all(pair_keys %in% colnames(df))) {
    df$.label <- if (length(pair_keys) == 1L) {
      as.character(df[[pair_keys]])
    } else {
      paste(df[[pair_keys[1]]], df[[pair_keys[2]]], sep = " - ")
    }
    top <- head(df[order(df$.p, na.last = TRUE), ], label_top)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = top,
                                        ggplot2::aes(label = .data$.label),
                                        size = 3, max.overlaps = Inf,
                                        show.legend = FALSE)
    } else {
      p <- p + ggplot2::geom_text(data = top,
                                  ggplot2::aes(label = .data$.label),
                                  size = 3, hjust = -0.1, vjust = -0.1,
                                  show.legend = FALSE)
    }
  }
  p
}
