library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(dplyr)
library(ggrastr)
library(pheatmap)
library(Matrix)
library(ggraph)
library(circlize)
library(scales)
library(ComplexHeatmap)
library(parallel)
library(RANN)
library(gridExtra)


# Functions

generate_sim = function(close_ratio = 0.7,  
                        n_types = 10, 
                        max_loc = 800,
                        n_perm = 100,
                        n_cells = 1500, 
                        test_type = "circle", # "line", "circle", "distribute"
                        distance_param = 50,  
                        seed=1234){
  if(test_type=="distribute"){
    
    set.seed(seed)  

    x_coords <- runif(n_cells, min = 0, max = max_loc)
    y_coords <- runif(n_cells, min = 0, max = max_loc)

    cell_types <- sample(paste0("cell_type_", 1:n_types), n_cells, replace = TRUE)

    idx_type_1 <- which(cell_types == "cell_type_1")
    idx_type_2 <- which(cell_types == "cell_type_2")

    n_pairs <- round(length(idx_type_2) * close_ratio)
    idx_close_1 <- sample(idx_type_1, n_pairs, replace = TRUE) 
    idx_close_2 <- sample(idx_type_2, n_pairs)

    angle_shift <- runif(n_pairs, 0, 2 * pi)
    x_coords[idx_close_2] <- x_coords[idx_close_1] + distance_param * cos(angle_shift) + rnorm(n_pairs, mean = 0, sd = distance_param/5)
    y_coords[idx_close_2] <- y_coords[idx_close_1] + distance_param * sin(angle_shift) + rnorm(n_pairs, mean = 0, sd = distance_param/5)

    df <- data.frame(x = x_coords, y = y_coords, cell_type = cell_types) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = paste0("cell_type_", 1:n_types)))
    
  } else if(test_type=="circle"){
    
    set.seed(seed)

    n1_total = n2_total = ceiling(n_cells/n_types)
    n_other <- n_cells-n1_total-n2_total

    close_ratio_1 = close_ratio_2 = close_ratio

    center_x <- ceiling(max_loc/2)
    center_y <- ceiling(max_loc/2)
    radius <- ceiling(center_x/3) 

    n1_close <- round(n1_total * close_ratio_1) 

    theta_1_outer <- seq(0, 2*pi, length.out = round((n1_close / 4)))
    x_1_outer <- center_x + radius * cos(theta_1_outer)# + rnorm(length(theta_1_outer), mean = distance_param, sd = ceiling(distance_param/10))
    y_1_outer <- center_y + radius * sin(theta_1_outer)# + rnorm(length(theta_1_outer), mean = distance_param, sd = ceiling(distance_param/10))
    
    r_1_inner <- sqrt(runif(round(n1_close - round((n1_close / 4))), 0, radius^2)) 
    theta_1_inner <- runif(round(n1_close - round((n1_close / 4))), 0, 2*pi)
    x_1_inner <- center_x + r_1_inner * cos(theta_1_inner)
    y_1_inner <- center_y + r_1_inner * sin(theta_1_inner)

    n1_random <- n1_total - length(x_1_inner) - length(x_1_outer) 
    x_1_random <- runif(n1_random, min = 0, max = max_loc)
    y_1_random <- runif(n1_random, min = 0, max = max_loc)

    x_1 <- c(x_1_inner, x_1_outer, x_1_random)
    y_1 <- c(y_1_inner, y_1_outer, y_1_random)

    n2_close <- round(n2_total * close_ratio_2) 

    theta_1_outer <- seq(0, 2*pi, length.out = round(n2_close))
    x_2_close <- center_x + (radius+distance_param) * cos(theta_1_outer) + rnorm(n2_close, mean = 0, sd = 5)
    y_2_close <- center_y + (radius+distance_param) * sin(theta_1_outer) + rnorm(n2_close, mean = 0, sd = 5)

    n2_random <- n2_total - n2_close 
    x_2_random <- runif(n2_random, min = 0, max = max_loc)
    y_2_random <- runif(n2_random, min = 0, max = max_loc)

    x_other <- runif(n_other, min = 0, max = max_loc)
    y_other <- runif(n_other, min = 0, max = max_loc)

    df <- data.frame(
      x = c(x_1, x_2_close, x_2_random, x_other),
      y = c(y_1, y_2_close, y_2_random, y_other),
      cell_type = c(rep("cell_type_1", n1_total), 
                    rep("cell_type_2", n2_close), 
                    rep("cell_type_2", n2_random), 
                    sample(paste0("cell_type_", 3:n_types), n_other, replace = TRUE))
    ) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = paste0("cell_type_", 1:n_types)))
    
  } else if (test_type=="line"){
    
    set.seed(seed)

    n_cells_per_layer <- ceiling(n_cells/n_types * max(close_ratio,0.01))
    n_layers <- 2 
    layer_spacing <- distance_param 
    x_start <- ceiling(max_loc/2)-ceiling(max_loc/4)
    y_start <- ceiling(max_loc/2)-50 
    x_end <- ceiling(max_loc/2)+ceiling(max_loc/4)
    step = (x_end-x_start)/n_cells_per_layer

    x_1 <- c()
    y_1 <- c()
    x_2 <- c()
    y_2 <- c()
    
    for (i in 0:(n_layers - 1)) {
      y_pos <- y_start + i * layer_spacing 
      
      if (i %% 2 == 0) {
        x_layer <- seq(x_start, x_start + n_cells_per_layer * step, length.out = n_cells_per_layer) + rnorm(n_cells_per_layer, mean = 0, sd = 1)
        y_layer <- rep(y_pos, n_cells_per_layer) + rnorm(n_cells_per_layer, 0, 2)

        n_close <- round(n_cells_per_layer * close_ratio)
        x_1 <- c(x_1, x_layer)
        y_1 <- c(y_1, y_layer)
      } else {
        x_layer <- seq(x_start, x_start + n_cells_per_layer * step, length.out = n_cells_per_layer) + rnorm(n_cells_per_layer, mean = 0, sd = 1)
        y_layer <- rep(y_pos, n_cells_per_layer) + rnorm(n_cells_per_layer, 0, 2)

        n_close <- round(n_cells_per_layer * close_ratio)
        idx_close <- sample(1:n_cells_per_layer, n_close, replace = TRUE)
        angle_shift <- runif(n_close, 0, 2 * pi)
        x_layer[idx_close] <- x_1[idx_close]
        y_layer[idx_close] <- y_1[idx_close] + distance_param + rnorm(n_close, mean = 0, sd = 1)
        
        x_2 <- c(x_2, x_layer)
        y_2 <- c(y_2, y_layer)
      }
    }

    n1_random = n2_random = ceiling(ceiling(n_cells/n_types) * (1-close_ratio))
    x_1_random <- runif(n1_random, min = 0, max = max_loc)
    y_1_random <- runif(n1_random, min = 0, max = max_loc)
    x_2_random <- runif(n2_random, min = 0, max = max_loc)
    y_2_random <- runif(n2_random, min = 0, max = max_loc)
    n_other <- n_cells - n1_random - n2_random
    x_other <- runif(n_other, min = 0, max = max_loc)
    y_other <- runif(n_other, min = 0, max = max_loc)

    df <- data.frame(
      x = c(x_1, x_1_random, x_2, x_2_random, x_other),
      y = c(y_1, y_1_random, y_2, y_2_random, y_other),
      cell_type = c(rep("cell_type_1", sum(length(x_1),length(x_1_random))), 
                    rep("cell_type_2", sum(length(x_2),length(x_2_random))), 
                    sample(paste0("cell_type_", 3:n_types), n_other, replace = TRUE))
    ) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = paste0("cell_type_", 1:n_types)))
    
  }
  return(df)
}


compute_count <- function(adj, int_clust_row, int_clust_col, n_cls, cluster_data, transformation = TRUE) {

  counts <- matrix(0, nrow = n_cls, ncol = n_cls, 
                   dimnames = list(paste0("Cluster", levels(cluster_data)), paste0("Cluster", levels(cluster_data))))
  
  for (i in paste0("Cluster", levels(cluster_data))) {
    cluster_rows <- which(int_clust_row == i)
    
    for (j in paste0("Cluster", levels(cluster_data))) {
      cluster_cols <- which(int_clust_col == j)
      neighbors <- adj[cluster_rows, cluster_cols, drop = FALSE]
      if(transformation){
        counts[i, j] <- sum(neighbors)
      } else {
        counts[i, j] <- sum(neighbors == 1) 
      }
    }
  }
  return(counts)
}

permute_clusters <- function(adj, int_clust, n_cls, cluster_data, transformation) {
  int_clust_row <- sample(int_clust) 
  int_clust_col <- sample(int_clust) 
  compute_count(adj, int_clust_row = int_clust_row, int_clust_col = int_clust_col, n_cls, cluster_data, transformation) 
}

calc_co_occurrence_for_radius <- function(
    seurat_obj, 
    radius, 
    sample_key,
    cluster_key,
    k = 30
) {

  all_clusters <- levels(factor(seurat_obj@meta.data[[cluster_key]]))

  co_occur_count <- matrix(
    0,
    nrow = length(all_clusters),
    ncol = length(all_clusters),
    dimnames = list(
      paste0("Cluster", all_clusters),
      paste0("Cluster", all_clusters)
    )
  )

  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]])) %>%
      tibble::column_to_rownames(var = "cell") %>%
      dplyr::mutate(
        cluster = seurat_obj@meta.data[
          seurat_obj@meta.data[, sample_key] == name,
          cluster_key
        ]
      )

    res <- RANN::nn2(
      data       = coords[, 1:2],   
      query      = coords[, 1:2],   
      searchtype = "radius",
      radius     = radius,
      k          = k
    )
    
    clusters_vec <- coords$cluster

    for (i in seq_len(nrow(coords))) {
      neighbors_i <- res$nn.idx[i, ]

      neighbors_i <- neighbors_i[neighbors_i != i & neighbors_i > 0]

      c_i <- paste0("Cluster", clusters_vec[i])
      
      if (length(neighbors_i) > 0) {
        c_neighbors <- paste0("Cluster", clusters_vec[neighbors_i])
        tab <- table(c_neighbors)
        for (cn in names(tab)) {
          co_occur_count[c_i, cn] <- co_occur_count[c_i, cn] + tab[[cn]]
        }
      }
    }
  }

  ratio_mat <- compute_co_occurrence_ratio(co_occur_count)

  return(list(
    co_occur_count = co_occur_count,
    ratio_mat      = ratio_mat
  ))
}


compute_co_occurrence_ratio <- function(co_occur_count) {
  rn <- rownames(co_occur_count)
  cn <- colnames(co_occur_count)
  row_sums <- rowSums(co_occur_count)
  col_sums <- colSums(co_occur_count)
  total_sum <- sum(co_occur_count)
  
  ratio_mat <- matrix(
    0, 
    nrow = nrow(co_occur_count), 
    ncol = ncol(co_occur_count),
    dimnames = list(rn, cn)
  )
  
  for (i in seq_len(nrow(co_occur_count))) {
    for (j in seq_len(ncol(co_occur_count))) {
      p_exp_cond <- co_occur_count[i, j] / row_sums[i]   # p(exp=j | cond=i)
      p_exp      <- col_sums[j]       / total_sum        # p(exp=j)
      ratio_mat[i, j] <- p_exp_cond / p_exp
    }
  }
  return(ratio_mat)
}

search_interaction_spot = function(seurat_object,
                                   fov,
                                   radius,
                                   n_min,
                                   neighbors.k = 200,
                                   cell_id = cell_id, # anchor cells
                                   cluster_col = cluster_col,
                                   target_cluster = target_cluster # cell clusters assumed to be interacted
){
  
  coords <- seurat_object[[fov]]$centroids@coords %>%
    as.data.frame() %>%
    dplyr::mutate(cell = Cells(seurat_object[[fov]])) %>%
    dplyr::filter(cell %in% cell_id)
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[, c("x", "y")])
  dim(coords)

  res_nn2 <- RANN::nn2(
    data       = coords,      
    query      = coords,      
    searchtype = "radius",
    radius     = radius,
    k          = neighbors.k
  )
  # (1) Create a ‘cell ID column’ in coords (use the row name if there is one)
  coords_df_ <- coords %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "cell_id")
  
  # (2) Create a graph based on nn.idx that shows ‘cells close to each other’
  # Draw an edge between the cell in row i and the cells listed in nn.idx[i, ]
  # (excluding 0)
  edges <- vector("list", length = nrow(res_nn2$nn.idx))  # 一旦リストで準備
  num_clus <- vector("list", length = nrow(res_nn2$nn.idx))
  for (i in seq_len(nrow(res_nn2$nn.idx))) {
    # Nearest index (excluding 0)
    neighbors_i <- res_nn2$nn.idx[i, res_nn2$nn.idx[i, ] != 0]
    # Add an edge from i to neighbors_i
    # Use c(rbind(...)) or lapply to make multiple pairs of (i, neighbors_i)
    edges[[i]] <- c(rbind(i, neighbors_i))
  }

  edges_vec <- unlist(edges, use.names = FALSE)
  
  # (3) Creating an igraph object
  # Assuming an undirected graph (directed=FALSE)
  g <- igraph::graph(edges = edges_vec, directed = FALSE)
  
  # (4) Search for a connected component to obtain a cluster ID
  comp <- igraph::components(g, mode = "weak")
  # component numbers (1, 2, 3, ...) to which each vertex (cell) belongs is stored in comp$membership.
  
  cluster_id <- comp$membership
  
  # (5) Add a cluster_id column to coords_df
  coords_df_$cluster_id <- cluster_id

  coords_df_ <- coords_df_ %>%
    dplyr::left_join(.,seurat_object@meta.data[,c("cell",cluster_col)] %>% dplyr::rename(cell_id = cell),by="cell_id") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(x_min = min(x),
                  x_max = max(x),
                  y_min = min(y),
                  y_max = max(y),
                  n_all_cells = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!!sym(cluster_col) %in% target_cluster) %>%
    dplyr::filter(n_all_cells > n_min) %>%
    dplyr::mutate(sample_id = fov,
                  cluster_id = paste0(fov,":cluster",cluster_id))
  
  return(coords_df_)
}

nhood_enrichment.Seurat <- function(seurat_obj, 
                             cluster_key, 
                             neighbors.k = 30, 
                             connectivity_key = "nn", 
                             transformation = TRUE,
                             n_perms = 100, seed = 1938493, n_jobs = 4) {
  
  set.seed(seed)
  
  if (!cluster_key %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cluster key", cluster_key, "not found in meta.data"))
  }
  cluster_data <- seurat_obj@meta.data[[cluster_key]]
  int_clust <- paste0("Cluster", cluster_data)
  
  all_nn <- list()  
  all_snn <- list() 
  cell_id = vector()
  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]]))
    cells <- coords$cell
    rownames(coords) <- cells
    coords <- as.matrix(coords[, c("x", "y")])
    neighbors <- FindNeighbors(coords, k.param = neighbors.k, verbose = FALSE)
    all_nn[[name]] <- neighbors$nn
    all_snn[[name]] <- neighbors$snn
    cell_id = c(cell_id,cells)
  }

  if(connectivity_key == "nn") {
    adj <- bdiag(all_nn)
  } else {
    adj <- bdiag(all_snn)
  }
  rownames(adj) = colnames(adj) = cell_id
  
  # --- normalize ---
  if (transformation) {
    degrees <- Matrix::colSums(adj) + 1 
    adj <- adj / degrees  
  }
  
  cluster_data = factor(cluster_data)
  n_cls <- length(levels(cluster_data))

  count <- compute_count(adj, int_clust_row = int_clust, int_clust_col = int_clust, n_cls, cluster_data, transformation)
  
  perform_permutations <- function(adj, int_clust, n_cls, cluster_data, n_perms, n_jobs, transformation) {
    
    cl <- makeCluster(n_jobs)
    on.exit(stopCluster(cl))  
    
    clusterExport(
      cl,
      c("compute_count", "permute_clusters", "adj", "int_clust", "n_cls", "cluster_data", "transformation"),
      envir = environment()
    )
    
    while(TRUE) {
      tryCatch({
        perms <- parLapply(cl, seq_len(n_perms), function(x) permute_clusters(adj, int_clust, n_cls, cluster_data, transformation))
        return(perms)
      }, error = function(e) {
        #message("Error occurred, retrying...")
        Sys.sleep(5) 
      })
    }
  }
  
  perms <- perform_permutations(adj, int_clust, n_cls, cluster_data, n_perms, n_jobs, transformation)
  
  compute_zscore <- function(counts, perms) {
    n_clusters <- ncol(counts)  
    zscore <- matrix(NA, nrow = n_clusters, ncol = n_clusters, 
                     dimnames = list(paste("Cluster", 1:n_clusters), paste("Cluster", 1:n_clusters)))
    
    for (i in seq_len(n_clusters)) {
      for (j in seq_len(n_clusters)) {
        perm_values <- sapply(perms, function(perm) perm[i, j]) 
        perm_mean <- mean(perm_values)
        perm_sd <- sd(perm_values)

        zscore[i, j] <- (counts[i, j] - perm_mean) / perm_sd
      }
    }
    
    return(zscore)
  }
  zscore <- compute_zscore(count, perms)
  
  perms_matrix <- do.call(rbind, perms)
  perm_means <- colMeans(perms_matrix)
  perm_sds <- apply(perms_matrix, 2, sd)
  zscore <- (count - perm_means) / perm_sds

  seurat_obj@misc[[paste0(cluster_key, "_nhood_enrichment")]] <- list(
    zscore = zscore,
    count = count
  )
  
  return(seurat_obj)
}

cooccur_local.Seurat = function(seurat_obj,
                                cluster_x,
                                cluster_y,
                                connectivity_key = "nn",
                                cluster_key = "seurat_clusters",
                                sample_key = "sample_id",
                                neighbors.k = 20, 
                                radius = 30,
                                maxnsteps= 15
){
  all_nn <- list() 
  all_snn <- list()  
  cell_id = vector()

  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]]))
    cells <- coords$cell
    rownames(coords) <- cells
    coords <- as.matrix(coords[, c("x", "y")])
    
    neighbors <- FindNeighbors(coords, k.param = neighbors.k, verbose = FALSE)
    
    all_nn[[name]] <- neighbors$nn
    all_snn[[name]] <- neighbors$snn
    cell_id = c(cell_id,cells)
  }
  
  if(connectivity_key == "nn") {
    adj <- bdiag(all_nn)
  } else {
    adj <- bdiag(all_snn)
  }
  rownames(adj) = colnames(adj) = cell_id
  
  local_score = vector()
  
  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]])) %>%
      tibble::column_to_rownames(var = "cell") %>%
      dplyr::mutate(
        cluster = seurat_obj@meta.data[
          seurat_obj@meta.data[, sample_key] == name,
          cluster_key
        ]
      )

    res_nn2 <- RANN::nn2(
      data       = coords[,1:2],      
      query      = coords[,1:2],      
      searchtype = "radius",
      radius     = radius,
      k          = neighbors.k    
    )

    
    local_score_ <- numeric(nrow(coords))  
    names(local_score_) <- rownames(coords)
    
    for (i in seq_len(nrow(coords))) {
      neighbors_i <- res_nn2$nn.idx[i, ]
      neighbors_i <- neighbors_i[neighbors_i > 0]
      neighbors_i <- neighbors_i[neighbors_i != i]
      
      if (length(neighbors_i) == 0) {
        local_score_[i] <- 0
        next
      }
      
      cluster_vec <- coords$cluster
      c_vec <- as.character(cluster_vec[neighbors_i])
      
      cond_present <- any(c_vec == cluster_y) & any(c_vec == cluster_x)
      
      if (cond_present) {
        local_score_[i] <- 1
      } else {
        local_score_[i] <- 0
      }
    }
    local_score = c(local_score, local_score_)
  }
  
  diffuse_step <- function(adj, local_score) {
    a <- adj  
    degrees <- Matrix::colSums(a) + 1 
    s_norm <- matrix(local_score) / degrees 
    res <- (a %*% s_norm) + s_norm  
    return(as.matrix(res))
  }  
  
  nsteps=NULL
  for (i in seq_len(maxnsteps)) { 
    s <- diffuse_step(adj, 
                      local_score)
    medkurt <- moments::kurtosis(prop.table(s, 2))
    if (is.null(nsteps)) { 
      prevmedkurt <- medkurt
      if(is.nan(medkurt)) {
        break
      }
      if (prevmedkurt - medkurt < 3 & i > 3) { 
        message(glue::glue('stopping after {i} steps'))
        break 
      }            
    } else if (i == nsteps) {
      break
    }
  }  
  
  df = data.frame(cooccur_local_scores = as.numeric(s)) %>%
    magrittr::set_colnames(paste0("cooccur_local_", cluster_x, "_", cluster_y))
  rownames(df) = names(local_score)
  return(df)
}

cooccur_local = function(df,
                         cluster_x,
                         cluster_y,
                         connectivity_key = "nn",
                         neighbors.k = 20, 
                         radius = 30,
                         maxnsteps= 1 
){
  if(class(df)=="Seurat"){
    return(cooccur_local.Seurat(df, cluster_x, cluster_y, connectivity_key, neighbors.k, radius, maxnsteps))
  }
  coords <- df %>%
    dplyr::mutate(cell = rownames(.))
  cell_id <- coords$cell
  rownames(coords) = coords$cell
  coords <- as.matrix(coords[, c("x", "y")])
  
  # transformation: advised for analytic p-value calculation.

  neighbors <- FindNeighbors(coords, k.param = neighbors.k, verbose = FALSE)
  
  if(connectivity_key == "nn") {
    adj <- neighbors$nn
  } else {
    adj <- neighbors$snn
  }
  rownames(adj) = colnames(adj) = cell_id

  coords <- df %>%
    dplyr::mutate(cell = rownames(.)) %>%
    dplyr::rename( cluster = cell_type)

  res_nn2 <- RANN::nn2(
    data       = coords[,1:2],  
    query      = coords[,1:2], 
    searchtype = "radius",
    radius     = radius,
    k          = neighbors.k    
  )

  local_score_ <- numeric(nrow(coords))
  names(local_score_) <- rownames(coords)
  
  for (i in seq_len(nrow(coords))) {
    neighbors_i <- res_nn2$nn.idx[i, ]
    neighbors_i <- neighbors_i[neighbors_i > 0]
    neighbors_i <- neighbors_i[neighbors_i != i]
    
    if (length(neighbors_i) == 0) {
      local_score_[i] <- 0
      next
    }

    cluster_vec <- coords$cluster
    c_vec <- as.character(cluster_vec[neighbors_i])
    cond_present <- any(c_vec == cluster_y) & any(c_vec == cluster_x)

    if (cond_present) {
      local_score_[i] <- 1
    } else {
      local_score_[i] <- 0
    }
  }
  local_score = local_score_
  
  if(maxnsteps == 0) {
    df_ = data.frame(cooccur_local_scores = as.numeric(local_score)) %>%
      magrittr::set_colnames(paste0("cooccur_local_", cluster_x, "_", cluster_y))
    rownames(df_) = names(local_score)
    return(df_)
  } else {
    diffuse_step <- function(adj, local_score) {
      a <- adj  
      degrees <- Matrix::colSums(a) + 1 
      s_norm <- matrix(local_score) / degrees  
      res <- (a %*% s_norm) + s_norm  
      return(as.matrix(res)) 
    }  
    
    nsteps=NULL
    for (i in seq_len(maxnsteps)) { 
      s <- diffuse_step(adj, 
                        local_score)
      medkurt <- moments::kurtosis(prop.table(s, 2))
      if (is.null(nsteps)) { 
        prevmedkurt <- medkurt
        if(is.nan(medkurt)) {
          break
        }
        if (prevmedkurt - medkurt < 3 & i > 3) { 
          #message(glue::glue('stopping after {i} steps'))
          break 
        }            
      } else if (i == nsteps) {
        break
      }
    }  
    df_ = data.frame(cooccur_local_scores = as.numeric(s)) %>%
      magrittr::set_colnames(paste0("cooccur_local_", cluster_x, "_", cluster_y))
    rownames(df_) = names(local_score)
    return(df_)
  }
}



nhood_enrichment <- function(df, 
                             cluster_key, 
                             neighbors.k = 30, 
                             connectivity_key = "nn", 
                             transformation = TRUE,
                             n_perms = 100, seed = 1938493, n_jobs = 4) {
  if(class(df)=="Seurat"){
    return(nhood_enrichment.Seurat(df, cluster_key, neighbors.k, connectivity_key, transformation, n_perms, seed, n_jobs))
  }
  set.seed(seed)
  
  cluster_data <- df[,cluster_key]
  int_clust <- paste0("Cluster", cluster_data)

  coords <- df %>%
    dplyr::mutate(cell = rownames(.))
  cell_id <- coords$cell
  rownames(coords) = coords$cell
  coords <- as.matrix(coords[, c("x", "y")])

  neighbors <- FindNeighbors(coords, k.param = neighbors.k, verbose = FALSE)
  
  if(connectivity_key == "nn") {
    adj <- neighbors$nn
  } else {
    adj <- neighbors$snn
  }
  rownames(adj) = colnames(adj) = cell_id
  
  # --- normalize ---
  if (transformation) {
    degrees <- Matrix::colSums(adj) + 1 
    adj <- adj / degrees 
  }
  
  cluster_data = factor(cluster_data)
  n_cls <- length(levels(cluster_data))
  
  count <- compute_count(adj, int_clust_row = int_clust, int_clust_col = int_clust, n_cls, cluster_data, transformation)
  
  perform_permutations <- function(adj, int_clust, n_cls, cluster_data, n_perms, n_jobs, transformation) {
    
    cl <- makeCluster(n_jobs)
    on.exit(stopCluster(cl)) 

    clusterExport(
      cl,
      c("compute_count", "permute_clusters", "adj", "int_clust", "n_cls", "cluster_data", "transformation"),
      envir = environment()
    )

    while(TRUE) {
      tryCatch({
        perms <- parLapply(cl, seq_len(n_perms), function(x) permute_clusters(adj, int_clust, n_cls, cluster_data, transformation))
        return(perms)
      }, error = function(e) {
        #message("Error occurred, retrying...")
        Sys.sleep(5)
      })
    }
  }
  
  perms <- perform_permutations(adj, int_clust, n_cls, cluster_data, n_perms, n_jobs, transformation)
  
  compute_zscore <- function(counts, perms) {
    n_clusters <- ncol(counts) 
    zscore <- matrix(NA, nrow = n_clusters, ncol = n_clusters, 
                     dimnames = list(paste("Cluster", 1:n_clusters), paste("Cluster", 1:n_clusters)))
    
    for (i in seq_len(n_clusters)) {
      for (j in seq_len(n_clusters)) {
        perm_values <- sapply(perms, function(perm) perm[i, j]) 
        perm_mean <- mean(perm_values)
        perm_sd <- sd(perm_values)
        zscore[i, j] <- (counts[i, j] - perm_mean) / perm_sd
      }
    }
    
    return(zscore)
  }
  zscore <- compute_zscore(count, perms)
  
  perms_matrix <- do.call(rbind, perms)
  perm_means <- colMeans(perms_matrix)
  perm_sds <- apply(perms_matrix, 2, sd)
  zscore <- (count - perm_means) / perm_sds
  
  return(list(
    zscore = zscore,
    count = count
  ))
  
}

manual_colors = c(
  "0" = "#E41A1C",
  "1" = "#377EB8",
  "2" = "#4DAF4A",
  "3" = "#984EA3",
  "4" = "#FF7F00",
  "5" = "#FFFF33",
  "6" = "#A65628",
  "7" = "lightgrey",
  "8" = "#999999",
  "9" = "#66C2A5",
  "10" = "#67000D",
  "11" = "#8DA0CB",
  "12" = "#FFD92F",
  "13" = "#A6D854",
  "14" = "#E78AC3",
  "15" = "#FC8D62",
  "16" = "darkgrey",
  "17" = "#FEB24C",
  "18" = "#377EB8",
  "19" = "lightblue",
  "20" = "#FDE0EF",
  "21" = "#B8E186",
  "22" = "#66C2A5",
  "23" = "#A6D855",
  "24" = "#E78AC4",
  "25" = "#FC8D63",
  "26" = "brown",
  "27" = "#FEB25C",
  "28" = "#377EB9",
  "29" = "lightgreen",
  "30" = "#FDE1EF"
)

