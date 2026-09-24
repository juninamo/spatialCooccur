# ---- Unsupervised, transcript-level co-localization modules ---------------
#
# pcf_matrix() needs gene sets chosen in advance (e.g. cell-type markers),
# and its results are only as good as those sets. The functions below work
# gene by gene: colocalization_gene_matrix() measures, for every pair of
# genes, how many transcript pairs lie within `radius` of each other compared
# with label shuffling over fixed transcript positions (the transcript
# analogue of log2 O/E); colocalization_modules() clusters genes that
# co-localize with one another; module_enrichment() interprets the modules
# with any list of gene sets (pathways, cell-type markers).

#' Gene-by-gene co-localization of transcripts
#'
#' **Experimental.** For genes \eqn{a, b}, counts transcript pairs within
#' distance `radius`,
#' \deqn{P_{ab} = \sum_u Y_a(u)\,(K_r * Y_b)(u),}
#' where \eqn{Y} are binned counts and \eqn{K_r} is a disc of radius `radius`
#' (FFT convolution, one per gene), and compares it with its expectation
#' when gene labels are shuffled over the fixed transcript positions,
#' \deqn{E_{ab} = \frac{n_a n_b}{N(N-1)} P_{\mathrm{all}},}
#' with \eqn{n} the transcript counts, \eqn{N} their total and
#' \eqn{P_{\mathrm{all}}} the number of all transcript pairs within the
#' radius (self pairs excluded). `log2((P + c) / (E + c))` is 0 without
#' co-localization; cellularity and tissue shape cancel, as in the relative
#' pair correlation.
#'
#' @param binned Output of [bin_transcripts()].
#' @param genes Genes to use (default: all genes with at least `min_count`
#'   transcripts).
#' @param radius Distance in um within which transcripts count as a pair.
#' @param min_count Minimum number of transcripts per gene.
#' @param top_n Keep at most this many genes (the most abundant).
#' @param pseudocount Added to observed and expected pair counts.
#'
#' @return A symmetric gene x gene matrix of log2 O/E with attributes
#'   `n_transcripts` (per gene) and `radius`.
#' @export
#' @examples
#' tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
#' b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
#' M <- colocalization_gene_matrix(b, radius = 12)
#' round(M[1:4, 1:4], 2)
colocalization_gene_matrix <- function(binned, genes = NULL, radius = 20, min_count = 50,
                                       top_n = NULL, pseudocount = 1) {
  if (!inherits(binned, "binned_transcripts")) stop("`binned` must come from bin_transcripts().")
  g <- binned$grid; nx <- g$nx; ny <- g$ny; d <- g$bin_size
  cnt <- binned$counts
  tot <- Matrix::colSums(cnt)
  if (is.null(genes)) genes <- colnames(cnt)[tot >= min_count]
  genes <- intersect(genes, colnames(cnt)[tot > 0])
  if (!is.null(top_n) && length(genes) > top_n) genes <- genes[order(-tot[genes])][seq_len(top_n)]
  if (length(genes) < 2) stop("fewer than two genes pass the filters.")
  inside <- as.numeric(binned$coords$in_tissue)
  Y <- cnt[, genes, drop = FALSE] * inside
  # disc kernel on the bin grid, embedded in a zero-padded FFT grid
  rb <- ceiling(radius / d)
  px <- stats::nextn(nx + rb); py <- stats::nextn(ny + rb)
  lag <- function(p) { k <- 0:(p - 1); ifelse(k < p / 2, k, k - p) }
  kern <- (d^2 * outer(lag(px)^2, lag(py)^2, "+")) <= radius^2
  fk <- stats::fft(kern * 1)
  conv <- function(v) {
    m <- matrix(0, px, py); m[seq_len(nx), seq_len(ny)] <- matrix(v, nx, ny)
    Re(stats::fft(stats::fft(m) * fk, inverse = TRUE))[seq_len(nx), seq_len(ny)] / (px * py)
  }
  G <- length(genes)
  P <- matrix(0, G, G, dimnames = list(genes, genes))
  for (j in seq_len(G)) {
    sj <- as.numeric(conv(as.numeric(Y[, j])))
    P[, j] <- as.numeric(Matrix::crossprod(Y, sj))
  }
  P <- (P + t(P)) / 2
  n <- Matrix::colSums(Y)
  diag(P) <- diag(P) - n                                # remove each transcript paired with itself
  yall <- Matrix::rowSums(cnt * inside)
  N <- sum(yall)
  Pall <- sum(yall * as.numeric(conv(yall))) - N
  E <- outer(n, n) * Pall / (N * (N - 1))
  diag(E) <- n * (n - 1) * Pall / (N * (N - 1))
  M <- log2((P + pseudocount) / (E + pseudocount))
  attr(M, "n_transcripts") <- n
  attr(M, "radius") <- radius
  M
}

#' Cluster genes into co-localization modules
#'
#' **Experimental.** Hierarchical clustering of the gene x gene log2 O/E
#' matrix from [colocalization_gene_matrix()]: genes whose transcripts lie
#' near one another (high mutual O/E) end up in the same module. The
#' distance between genes is `max(M) - M` with average linkage.
#'
#' @param M Output of [colocalization_gene_matrix()].
#' @param n_modules Number of modules (`cutree(k = )`). If `NULL`, the tree
#'   is cut at height `max(M) - min_oe`, i.e. genes join a module when their
#'   average mutual log2 O/E exceeds `min_oe`.
#' @param min_oe Threshold used when `n_modules` is `NULL`.
#' @param min_size Modules with fewer genes are labelled `NA`.
#'
#' @return A list with `modules` (data.frame: gene, module, connectivity =
#'   mean log2 O/E with the other genes of its module), `summary`
#'   (data.frame: module, size, mean_oe, top_genes) and the `tree`.
#' @export
#' @examples
#' tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
#' b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
#' mods <- colocalization_modules(colocalization_gene_matrix(b, radius = 12), n_modules = 3)
#' mods$summary
colocalization_modules <- function(M, n_modules = NULL, min_oe = 0.5, min_size = 3) {
  D <- stats::as.dist(max(M) - M)
  tr <- stats::hclust(D, method = "average")
  cl <- if (!is.null(n_modules)) stats::cutree(tr, k = n_modules) else stats::cutree(tr, h = max(M) - min_oe)
  sz <- table(cl)
  cl[cl %in% as.integer(names(sz)[sz < min_size])] <- NA
  # renumber by size
  lv <- names(sort(table(cl), decreasing = TRUE))
  mod <- ifelse(is.na(cl), NA, paste0("M", match(as.character(cl), lv)))
  genes <- rownames(M)
  conn <- vapply(seq_along(genes), function(i) {
    if (is.na(mod[i])) return(NA_real_)
    o <- setdiff(which(mod == mod[i]), i); if (!length(o)) NA_real_ else mean(M[i, o])
  }, 0)
  modules <- data.frame(gene = genes, module = mod, connectivity = conn, stringsAsFactors = FALSE)
  modules <- modules[order(modules$module, -modules$connectivity), ]
  summ <- do.call(rbind, lapply(split(modules, modules$module), function(m) {
    idx <- match(m$gene, genes); sub <- M[idx, idx, drop = FALSE]
    data.frame(module = m$module[1], size = nrow(m), mean_oe = mean(sub[upper.tri(sub)]),
               top_genes = paste(utils::head(m$gene, 8), collapse = ", "))
  }))
  summ <- summ[order(as.integer(sub("M", "", summ$module))), ]
  rownames(summ) <- NULL
  list(modules = modules, summary = summ, tree = tr)
}

#' Interpret co-localization modules with gene sets
#'
#' **Experimental.** Over-representation of each gene set (pathways,
#' cell-type markers, ...) in each module by the hypergeometric test,
#' against the genes that were clustered (`universe`), with
#' Benjamini-Hochberg correction over all module x set tests. Gene sets can
#' come from anywhere, e.g. `msigdbr::msigdbr()` or enrichR libraries.
#'
#' @param modules Output of [colocalization_modules()] (or its `modules`
#'   data.frame).
#' @param gene_sets Named list of character vectors.
#' @param universe Background genes (default: all clustered genes).
#' @param min_overlap Minimum overlap to report.
#'
#' @return A data.frame: module, gene_set, overlap, module_size, set_size
#'   (within the universe), odds_ratio, p, padj, genes.
#' @export
#' @examples
#' tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
#' b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
#' mods <- colocalization_modules(colocalization_gene_matrix(b, radius = 12), n_modules = 3)
#' truth <- attr(tx, "truth")
#' module_enrichment(mods, split(truth$genes, truth$set_of))
module_enrichment <- function(modules, gene_sets, universe = NULL, min_overlap = 1) {
  m <- if (is.data.frame(modules)) modules else modules$modules
  if (is.null(universe)) universe <- m$gene
  m <- m[!is.na(m$module) & m$gene %in% universe, ]
  U <- length(universe)
  gs <- lapply(gene_sets, function(s) intersect(unique(s), universe))
  gs <- gs[lengths(gs) > 0]
  res <- do.call(rbind, lapply(split(m$gene, m$module), function(mg) do.call(rbind, lapply(names(gs), function(k) {
    ov <- intersect(mg, gs[[k]]); a <- length(ov); K <- length(gs[[k]]); n <- length(mg)
    data.frame(gene_set = k, overlap = a, module_size = n, set_size = K,
               odds_ratio = (a + 0.5) * (U - K - n + a + 0.5) / ((n - a + 0.5) * (K - a + 0.5)),
               p = stats::phyper(a - 1, K, U - K, n, lower.tail = FALSE), genes = paste(ov, collapse = ", "))
  }))))
  res$module <- sub("\\.[0-9]+$", "", rownames(res)); rownames(res) <- NULL
  res$padj <- stats::p.adjust(res$p, method = "BH")
  res <- res[res$overlap >= min_overlap, c("module", "gene_set", "overlap", "module_size", "set_size", "odds_ratio", "p", "padj", "genes")]
  res[order(res$p), ]
}

#' Query enrichR for each co-localization module
#'
#' **Experimental.** Convenience wrapper around [enrichR::enrichr()] (needs
#' the enrichR package and an internet connection): submits the gene names
#' of every module to the Enrichr web service and returns the combined
#' results. Only gene names are sent. For offline use, or to test your own
#' pathways or cell-type markers, use [module_enrichment()].
#'
#' @param modules Output of [colocalization_modules()].
#' @param databases enrichR libraries, e.g. `"GO_Biological_Process_2023"`,
#'   `"CellMarker_2024"`, `"MSigDB_Hallmark_2020"`.
#' @param site Enrichr site passed to [enrichR::setEnrichrSite()].
#'
#' @return A data.frame with a `module` and a `database` column added to the
#'   enrichR output.
#' @export
module_enrichr <- function(modules, databases = c("GO_Biological_Process_2023", "CellMarker_2024"), site = "Enrichr") {
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("module_enrichr() requires the 'enrichR' package.")
  # enrichR sets its server options when it is attached
  if (is.null(getOption("enrichR.sites.base.address")) && !"package:enrichR" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("enrichR"))
  }
  suppressMessages(enrichR::setEnrichrSite(site))
  m <- if (is.data.frame(modules)) modules else modules$modules
  m <- m[!is.na(m$module), ]
  do.call(rbind, lapply(split(m$gene, m$module), function(g) {
    r <- enrichR::enrichr(g, databases)
    do.call(rbind, lapply(names(r), function(db) if (nrow(r[[db]])) cbind(module = m$module[match(g[1], m$gene)], database = db, r[[db]])))
  }))
}
