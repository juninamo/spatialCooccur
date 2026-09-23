# Tutorial for Spatial Neighborhood Analysis (SNA) & Spatial Co-localization Score (sCLS) Using Simulation Data

This article is a rendered copy of the Jupyter notebook
[`vignettes/SNA_tutorial_simulation.ipynb`](https://github.com/juninamo/spatialCooccur/blob/master/vignettes/SNA_tutorial_simulation.ipynb);
download it to run the code yourself.

**Author:**  
**Jun Inamo**  
*Computational Omics and Systems Immunology (COSI) Lab*  
*Division of Rheumatology and Center for Health AI*  
*University of Colorado School of Medicine, CO, USA*  
📧 <jun.inamo@cuanschutz.edu>

``` r

format(Sys.time(), '%d %B, %Y')
```

‘08 October, 2025’

## Spatial neighborhood analysis (SNA, cell type level analysis)

First generate dummy data where two cell types (cell_type_1 and
cell_type_2 in “cell_type” cik) are close to each other. Here we assume
that cell_type_1 and cell_type_2 are close to each other (concentric
cicle) by close_ratio. Make sure to assume \>10 cell types to see the
effect of neighborhood enrichment analysis. If this total cell types are
less than 10, the analysis may by susceptible to false positive results
because cell labels are shuffling in permutation test.

``` r

library(spatialCooccur)
library(patchwork)
library(ggplot2)
library(magrittr)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(ggrastr)
```

``` r

seed=1234
close_ratio=1  # proportion of cell_type_1 and cell_type_2 that are close to each other
n_types=15  # Number of cell types
max_loc=800 # Maximum x and y coordinates in space
n_cells=500  # Number of total cells
test_type="circle" # "line", "circle", "distribute"
distance_param=20  # distance between cell_type_1 and cell_type_2

df = generate_sim(close_ratio = close_ratio, 
                  n_types = n_types,  
                  max_loc = max_loc,
                  n_cells = n_cells,  
                  test_type = test_type,
                  distance_param = distance_param,  
                  seed=seed)
head(df) 
# x and y: coordinates
# cell_type: cell type
```

|     | x        | y        | cell_type   |
|-----|----------|----------|-------------|
|     | \<dbl\>  | \<dbl\>  | \<fct\>     |
| 1   | 355.4030 | 392.7360 | cell_type_1 |
| 2   | 490.8708 | 345.9954 | cell_type_1 |
| 3   | 451.1619 | 308.7717 | cell_type_1 |
| 4   | 501.4539 | 430.0083 | cell_type_1 |
| 5   | 280.3691 | 433.8681 | cell_type_1 |
| 6   | 389.7839 | 506.7382 | cell_type_1 |

A data.frame: 6 × 3 {.table .dataframe}

check the distribution of cell types in the space

``` r

cluster_colors = manual_colors
names(cluster_colors) = paste0("cell_type_",1:n_types)

g1 = ggplot(df, aes(x = x, y = y, color = cell_type)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = cluster_colors) +
  theme_void() +
  ggtitle("Simulated Cell Interaction Data")

g2 = ggplot(df %>% dplyr::mutate(cell_type = ifelse(cell_type %in% c("cell_type_1","cell_type_2"),as.character(cell_type),"others")), aes(x = x, y = y, fill = cell_type)) +
  geom_point(data = df %>% dplyr::mutate(cell_type = ifelse(cell_type %in% c("cell_type_1","cell_type_2"),as.character(cell_type),"others")) %>% dplyr::filter(cell_type %in% c("others")),
             alpha = 0.7, size = 2, shape = 21, stroke = 0.05, color = "black") +
  geom_point(data = df %>% dplyr::mutate(cell_type = ifelse(cell_type %in% c("cell_type_1","cell_type_2"),as.character(cell_type),"others")) %>% dplyr::filter(cell_type %in% c("cell_type_1","cell_type_2")),
             alpha = 0.8, size = 2, shape = 21, stroke = 0.05, color = "black") +
  labs(title = "Simulated Cell Interaction Data",
       subtitle = paste("n_types:", n_types, "| close_ratio:", close_ratio, "| max_loc:", max_loc, "| n_cells:", n_cells, "| distance_param:",distance_param)) +
  scale_color_manual(values = c("cell_type_1" = "red", 
                                "cell_type_2" = "blue",
                                "others" = "grey90")) +
  scale_fill_manual(values = c("cell_type_1" = "red", 
                               "cell_type_2" = "blue",
                               "others" = "grey90")) +
  theme_void()

  options(repr.plot.width=8, repr.plot.height=3)
g1|g2
```

``` output
“No shared levels found between `names(values)` of the manual scale and the
data's colour values.”
```

![](figures/SNA_tutorial_simulation/fig-01.png)

Run neighborhood enrichment analysis

``` r

n_perm = 100
neighbors.k_ = 30 # Number of neighbors to search

nhood_enrichment_res <- nhood_enrichment(
  df,
  cluster_key = "cell_type", 
  neighbors.k = neighbors.k_, 
  connectivity_key = "nn", 
  transformation = TRUE,
  n_perms = n_perm, seed = seed, n_jobs = 4
)
nhood_enrichment_res=nhood_enrichment_res$zscore

nhood_enrichment_res
```

|  | Clustercell_type_1 | Clustercell_type_2 | Clustercell_type_3 | Clustercell_type_4 | Clustercell_type_5 | Clustercell_type_6 | Clustercell_type_7 | Clustercell_type_8 | Clustercell_type_9 | Clustercell_type_10 | Clustercell_type_11 | Clustercell_type_12 | Clustercell_type_13 | Clustercell_type_14 | Clustercell_type_15 |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| Clustercell_type_1 | 17.9871780 | 4.7817679 | -3.65314566 | -2.12732491 | -3.19255339 | -2.29393557 | -3.5154081 | -2.78522419 | -3.38001771 | 1.301076 | -1.74187916 | -3.9894063 | -1.3912251 | -2.105634361 | -2.61864572 |
| Clustercell_type_2 | 7.6658961 | 10.8842376 | -2.69373489 | -0.39481708 | -1.36416923 | -2.98575264 | -2.2182959 | -3.27224251 | -2.12733430 | -1.017797 | -0.76764430 | -3.3963502 | -1.2794098 | -1.758017008 | -2.49274076 |
| Clustercell_type_3 | -2.9761619 | -1.0090318 | 1.55200614 | 1.54866861 | 1.82634482 | -0.04462895 | 1.3230457 | 0.10670713 | 0.03756225 | 2.303521 | 0.50802505 | -1.0687816 | -0.2602177 | -0.002965376 | 1.06972714 |
| Clustercell_type_4 | -1.6227139 | -0.1241113 | -0.18820728 | 2.22619881 | 1.09374148 | -1.21103296 | 0.3551040 | -0.66545726 | -0.38123679 | 1.966876 | 0.38450573 | -2.1929610 | -0.9272234 | 0.397422056 | -0.28458741 |
| Clustercell_type_5 | -3.4138188 | -1.5279760 | 0.11949763 | 1.15928393 | 3.43591067 | -0.23769514 | 1.3706482 | 0.30817463 | 0.73702581 | 1.805702 | -0.42712246 | -0.6147119 | 0.2317648 | 0.879774872 | 1.31679915 |
| Clustercell_type_6 | -0.6517928 | -1.3779308 | -0.74393169 | 0.36344031 | 0.91042137 | 1.10737342 | 0.3046584 | -0.30686414 | 0.35023999 | 2.314959 | -1.17527492 | -0.9991580 | -0.6887218 | -0.226417712 | 1.95624273 |
| Clustercell_type_7 | -3.3627457 | -1.3291273 | 0.30187928 | 0.85697197 | 1.77237156 | 0.31636230 | 3.3593937 | -0.43498178 | 0.56823891 | 2.060907 | -0.82550974 | -1.3659994 | -0.3794609 | 0.296483634 | 0.73622879 |
| Clustercell_type_8 | -2.6349217 | -2.6670134 | -0.16373606 | 0.76367570 | 1.78053052 | -0.55088908 | 0.3173527 | 1.79688893 | 0.03821519 | 2.226846 | 0.13466433 | -2.0951257 | 0.6769182 | 1.491857747 | 1.30491052 |
| Clustercell_type_9 | -2.7687376 | -1.3555109 | -0.70955906 | 0.25574694 | 1.26850775 | -0.01671934 | 1.1578975 | 0.03175541 | 2.14076649 | 2.033952 | -0.56776949 | -0.8820415 | 0.5456851 | -0.031767696 | 0.50685706 |
| Clustercell_type_10 | 0.5652055 | -2.0138750 | -1.16961727 | 0.15320208 | 0.87147238 | -0.21664598 | 0.5088938 | -0.35731577 | -0.88522898 | 3.068639 | -0.44068579 | -1.7898222 | -0.9956563 | 0.278669313 | 0.85695169 |
| Clustercell_type_11 | -0.5579479 | 0.4252563 | -0.01873493 | 0.91857621 | -0.03910266 | -2.16110188 | -0.1223153 | 0.38885975 | -0.85715327 | 2.379773 | 1.18954946 | -1.6673539 | -0.1508552 | 0.720147878 | 0.01940884 |
| Clustercell_type_12 | -2.4749385 | -1.2260379 | -0.69172407 | 0.28653298 | 1.06433717 | -0.67824403 | 0.1032652 | -0.93980525 | 0.25904476 | 1.686511 | -0.12302033 | 0.5669522 | 0.4234686 | 0.395742783 | 0.71856560 |
| Clustercell_type_13 | -1.4962686 | -1.2253934 | -0.15155209 | 0.37352104 | 1.19832098 | -1.51234488 | 0.4175916 | 0.65498550 | -0.02784879 | 1.589921 | -0.32823256 | -1.2999869 | 2.7341507 | 0.587693518 | 0.86209196 |
| Clustercell_type_14 | -1.4226028 | -0.6171552 | -0.85396030 | 0.06652673 | 1.62351940 | -1.79007462 | 0.1438312 | 0.21505198 | -1.26452053 | 1.764798 | -0.35795765 | -1.4471873 | -0.7671178 | 1.826489379 | -0.22236356 |
| Clustercell_type_15 | -2.6830525 | -1.8023621 | -1.02035397 | 0.55174468 | 1.77494723 | -0.00690145 | 0.4204384 | 0.57575913 | -0.33875210 | 2.532334 | 0.03504458 | -1.5008791 | 0.2291114 | 0.717072639 | 2.32544618 |

A matrix: 15 × 15 of type dbl {.table .dataframe}

How to read the Z-score matrix • Row: Reference cell type (cell_type_i)
• Column: Nearby cell type (cell_type_j) • The value of
nhood_enrichment\[i, j\] is the Z-score of ‘how close cell_type_i is to
cell_type_j’

Therefore, • nhood_enrichment\[‘cell_type_1’, ‘cell_type_2’\] • How much
cell_type_2 is gathered near cell_type_1 •
nhood_enrichment\[‘cell_type_2’, ‘cell_type_1’\] • How many cell_type_1s
are near cell_type_2?

Appropriate selection • When comparing a two-way relationship: check
both nhood_enrichment\[‘cell_type_1’, ‘cell_type_2’\] and
nhood_enrichment\[‘cell_type_2’, ‘cell_type_1’\] • When investigating
only one-way relationships: • Is there a high concentration of
cell_type_2 around cell_type_1? → nhood_enrichment\[‘cell_type_1’,
‘cell_type_2’\] • Is there a high concentration of cell_type_1 around
cell_type_2? → nhood_enrichment\[‘cell_type_2’, ‘cell_type_1’\]

``` r

colnames(nhood_enrichment_res) = gsub("^Cluster","",colnames(nhood_enrichment_res))
rownames(nhood_enrichment_res) = gsub("^Cluster","",rownames(nhood_enrichment_res))
pval_mat <- 1 - pnorm(nhood_enrichment_res)

fdr_vec <- p.adjust(as.vector(pval_mat), method = "BH")
fdr_mat <- matrix(fdr_vec, nrow=nrow(nhood_enrichment_res), ncol=ncol(nhood_enrichment_res),
                  dimnames = dimnames(nhood_enrichment_res))

sig_mat <- ifelse(fdr_mat < 0.05, "**", ifelse(fdr_mat > 0.05 & fdr_mat < 0.1, "*", ""))

common_names <- intersect(rownames(nhood_enrichment_res), colnames(nhood_enrichment_res))
for (nm in common_names) {
  nhood_enrichment_res[nm, nm] <- NA
  sig_mat[nm, nm] <- ""
}

heatmap <- Heatmap(nhood_enrichment_res,
                   name = "Z-score",
                   col = colorRamp2(c(-2, 0, 2), c("#0072B5FF", "white", "#BC3C29FF")), 
                   show_row_names = TRUE, 
                   show_column_names = TRUE,  
                   cluster_rows = TRUE,  
                   cluster_columns = TRUE,  
                   #show_column_dend = FALSE,
                   #show_row_dend = FALSE,
                   row_title = "",  
                   column_title = "Spatial Neigborhood Enrichment by cell types",
                   rect_gp = gpar(col = "black", lwd = 0.3),
                   na_col = "black",          # make sure cell with same cell types in the row and column be NA
                   
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if(sig_mat[i, j] == "**") {
                       grid.text("**", 
                                 x = x,
                                 y = y - 0.2 * height,  
                                 gp = gpar(fontsize = 15, col = "white", fontface = "bold"))
                     }
                     if(sig_mat[i, j] == "*") {
                       grid.text("*", 
                                 x = x,
                                 y = y - 0.2 * height,  
                                 gp = gpar(fontsize = 15, col = "white", fontface = "bold"))
                     }
                   }
)

options(repr.plot.width=6, repr.plot.height=6)
draw(heatmap, 
     merge_legend = TRUE,
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
```

![](figures/SNA_tutorial_simulation/fig-02.png)

Run the same analysis with different distance parameters

``` r

set.seed(seed)
random_seeds <- sample(1000:9999, 100)

accuracy_df_all = data.frame()

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

accuracy_df <- foreach(
  distance_param = c(5,10,20,30,40,50,75,100),
  .combine = rbind, 
  .packages = c("ggplot2", "dplyr", "Seurat", "doParallel", "spatialCooccur")
) %:%
  foreach(
    test_type = c("circle"), 
    .combine = rbind
  ) %:%
  foreach(
    seed_ = random_seeds,
    .combine = rbind
  ) %dopar% {
    
    set.seed(seed_) 
    print(grep(seed_, random_seeds))
    
    df = generate_sim(close_ratio = close_ratio, 
                      n_types = n_types,  
                      max_loc = max_loc,
                      n_cells = n_cells,  
                      test_type = test_type,
                      distance_param = distance_param,  
                      seed=1234)
    
    nhood_enrichment_res <- nhood_enrichment(
      df,
      cluster_key = "cell_type", 
      neighbors.k = neighbors.k_, 
      connectivity_key = "nn", 
      transformation = TRUE,
      n_perms = n_perm, seed = seed_, n_jobs = 4
    )
    
    nhood_enrichment_res <- nhood_enrichment_res$zscore
    diag(nhood_enrichment_res) <- 0
    colnames(nhood_enrichment_res) <- gsub("^Cluster", "", colnames(nhood_enrichment_res))
    rownames(nhood_enrichment_res) <- gsub("^Cluster", "", rownames(nhood_enrichment_res))
    
    pval_mat <- 1 - pnorm(nhood_enrichment_res)
    fdr_vec <- p.adjust(as.vector(pval_mat), method = "BH")
    fdr_mat <- matrix(fdr_vec, nrow = nrow(nhood_enrichment_res), ncol = ncol(nhood_enrichment_res), dimnames = dimnames(nhood_enrichment_res))
    
    sig_mat <- ifelse(fdr_mat < 0.05, "**", ifelse(fdr_mat > 0.05 & fdr_mat < 0.1, "*", ""))
    
    result <- data.frame(
      test_type = test_type,
      seed = seed_,
      distance_param = distance_param,
      zscore = nhood_enrichment_res["cell_type_1", "cell_type_2"],
      zscore_false = nhood_enrichment_res["cell_type_3", "cell_type_4"]
    )
    
    return(result)
  }

stopCluster(cl)

accuracy_df <- as.data.frame(accuracy_df)
head(accuracy_df)
accuracy_df_all = rbind(accuracy_df_all,accuracy_df)

accuracy_df_all$n_types = n_types
accuracy_df_all$neighbors.k_ = neighbors.k_
accuracy_df_all$close_ratio = close_ratio
accuracy_df_all$n_perm = n_perm
accuracy_df_all$general_seed = seed
accuracy_df_all$max_loc = max_loc
accuracy_df_all$n_cells = n_cells
```

|     | test_type | seed    | distance_param | zscore   | zscore_false |
|-----|-----------|---------|----------------|----------|--------------|
|     | \<chr\>   | \<int\> | \<dbl\>        | \<dbl\>  | \<dbl\>      |
| 1   | circle    | 8451    | 5              | 7.099714 | 1.685421     |
| 2   | circle    | 9015    | 5              | 7.010090 | 1.658070     |
| 3   | circle    | 8161    | 5              | 7.032457 | 1.620338     |
| 4   | circle    | 9085    | 5              | 6.899482 | 1.601554     |
| 5   | circle    | 8268    | 5              | 6.983298 | 1.655557     |
| 6   | circle    | 1622    | 5              | 7.135838 | 1.624043     |

A data.frame: 6 × 5 {.table .dataframe}

``` r

summary_df <- accuracy_df_all %>%
  dplyr::mutate(zscore_false = abs(zscore_false)) %>%
  tidyr::pivot_longer(cols = c(zscore, zscore_false), names_to = "zscore_type", values_to = "zscore") %>%
  dplyr::group_by(zscore_type, test_type, distance_param) %>%
  dplyr::summarise(
    median_zscore = median(zscore),
    lower_ci = quantile(zscore, 0.025),
    upper_ci = quantile(zscore, 0.975)
  ) %>%
  dplyr::ungroup()

options(repr.plot.width=6, repr.plot.height=6)
ggplot(summary_df, aes(x = distance_param, y = median_zscore, color = zscore_type)) +
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) + 
  geom_point(size = 2) + 
  facet_wrap(test_type ~ .) + 
  labs(
    x = "Distance Parameter",
    y = "Z-score",
    title = "Z-score Median and Confidence Interval by Distance Parameter",
    subtitle = paste("n_types:", n_types, "| neighbors.k_:", neighbors.k_, "| close_ratio:", close_ratio, "| max_loc:", max_loc, "| n_cells:", n_cells)
  ) +
  scale_x_log10(breaks = c(unique(summary_df$distance_param))) + 
  scale_color_manual(values = c("zscore" = "red", "zscore_false" = "grey40")) +  
  geom_hline(yintercept = 0, linetype = "dashed") +  
  #geom_hline(yintercept = 1.96, linetype = "dashed", color = "grey70") +  
  geom_hline(yintercept = abs(qnorm((0.05/2/n_types),F)), linetype = "dashed", color = "grey70") + 
  theme_classic() +
  theme(
    strip.text = element_text(size = 14, face = "bold"), 
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1)
  )
```

``` output
`summarise()` has grouped output by 'zscore_type', 'test_type'. You can
override using the `.groups` argument.
“Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
ℹ Please use `linewidth` instead.”
```

![](figures/SNA_tutorial_simulation/fig-03.png)

## Spatial co-localization score (sCLA, cell-cell level analysis)

Again we generate dummy data where two cell types (cell_type_1 and
cell_type_2 in “cell_type” cik) are close to each other. Here we assume
that cell_type_1 and cell_type_2 are close to each other (concentric
cicle) by close_ratio.

``` r

seed=1234
close_ratio=1  # proportion of cell_type_1 and cell_type_2 that are close to each other
n_types=10  # Number of cell types
max_loc=800 # Maximum x and y coordinates in space
n_cells=500  # Number of total cells
test_type="circle" # "line", "circle", "distribute"
distance_param = 20  # distance between cell_type_1 and cell_type_2

df = generate_sim(close_ratio = close_ratio, 
                  n_types = n_types,  
                  max_loc = max_loc,
                  n_cells = n_cells,  
                  test_type = test_type,
                  distance_param = distance_param,  
                  seed=seed)
head(df) 
# x and y: coordinates
# cell_type: cell type
```

|     | x        | y        | cell_type   |
|-----|----------|----------|-------------|
|     | \<dbl\>  | \<dbl\>  | \<fct\>     |
| 1   | 445.1298 | 397.7724 | cell_type_1 |
| 2   | 437.2731 | 301.0822 | cell_type_1 |
| 3   | 301.2229 | 365.6020 | cell_type_1 |
| 4   | 335.8962 | 315.8329 | cell_type_1 |
| 5   | 352.9081 | 515.0693 | cell_type_1 |
| 6   | 322.6797 | 325.7101 | cell_type_1 |

A data.frame: 6 × 3 {.table .dataframe}

check the distribution of cell types in the space

``` r

cluster_colors = manual_colors
names(cluster_colors) = paste0("cell_type_",1:n_types)

g1 = ggplot(df, aes(x = x, y = y, color = cell_type)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = cluster_colors) +
  theme_void() +
  ggtitle("Simulated Cell Interaction Data")

g2 = ggplot(df %>% dplyr::mutate(cell_type = ifelse(cell_type %in% c("cell_type_1","cell_type_2"),as.character(cell_type),"others")), aes(x = x, y = y, fill = cell_type)) +
  geom_point(data = df %>% dplyr::mutate(cell_type = ifelse(cell_type %in% c("cell_type_1","cell_type_2"),as.character(cell_type),"others")) %>% dplyr::filter(cell_type %in% c("others")),
             alpha = 0.7, size = 2, shape = 21, stroke = 0.05, color = "black") +
  geom_point(data = df %>% dplyr::mutate(cell_type = ifelse(cell_type %in% c("cell_type_1","cell_type_2"),as.character(cell_type),"others")) %>% dplyr::filter(cell_type %in% c("cell_type_1","cell_type_2")),
             alpha = 0.8, size = 2, shape = 21, stroke = 0.05, color = "black") +
  labs(title = "Simulated Cell Interaction Data",
       subtitle = paste("n_types:", n_types, "| close_ratio:", close_ratio, "| max_loc:", max_loc, "| n_cells:", n_cells, "| distance_param:",distance_param)) +
  scale_color_manual(values = c("cell_type_1" = "red", 
                                "cell_type_2" = "blue",
                                "others" = "grey90")) +
  scale_fill_manual(values = c("cell_type_1" = "red", 
                               "cell_type_2" = "blue",
                               "others" = "grey90")) +
  theme_void()

options(repr.plot.width=8, repr.plot.height=4)
g1|g2
```

``` output
“No shared levels found between `names(values)` of the manual scale and the
data's colour values.”
```

![](figures/SNA_tutorial_simulation/fig-04.png)

Run co-localization analysis

``` r

cluster_x <- "cell_type_1"
cluster_y <- "cell_type_2"

radius_ = 30 # Radius to search for neighbors (µm, cell_type_2) around anchor cells (cell_type_1)
neighbors.k_ = 30 # Number of neighbors to search

cooccur_local_df <- cooccur_local(
  df,
  cluster_x        = cluster_x,
  cluster_y        = cluster_y,
  connectivity_key = "nn",
  neighbors.k      = neighbors.k_, 
  radius           = radius_,
  maxnsteps        = 15
)
summary(cooccur_local_df)
```

``` output
 cooccur_local_cell_type_1_cell_type_2
 Min.   :0.0000                       
 1st Qu.:0.0000                       
 Median :0.0000                       
 Mean   :0.1020                       
 3rd Qu.:0.2141                       
 Max.   :0.3982                       
```

Check how the co-localization score is distributed in the space

``` r

g3 = ggplot() +
  geom_point(df %>%
               dplyr::mutate(score = cooccur_local_df$cooccur_local_cell_type_1_cell_type_2) %>%
               dplyr::filter(!(cell_type %in% c("cell_type_1","cell_type_2"))), 
             mapping = aes(x = x, y = y, fill = score),
             alpha = 0.7, size = 2, shape = 21, stroke = 0.05, color = "black") +
  geom_point(df %>%
               dplyr::mutate(score = cooccur_local_df$cooccur_local_cell_type_1_cell_type_2) %>%
               dplyr::filter(cell_type %in% c("cell_type_1","cell_type_2")), 
             mapping = aes(x = x, y = y, fill = score),
             alpha = 0.7, size = 2, shape = 21, stroke = 0.05, color = "black") +
  theme_void() +
  scale_fill_gradient(low = "grey90", high = "red") +
  ggtitle("Simulated Cell Interaction Data")

options(repr.plot.width=12, repr.plot.height=4)
g1 | g2 | g3
```

``` output
“No shared levels found between `names(values)` of the manual scale and the
data's colour values.”
```

![](figures/SNA_tutorial_simulation/fig-05.png)

Check the relationship between the co-localization score and the
distance between cell_type_1 and cell_type_2

``` r

one_cells <- df %>% dplyr::filter(cell_type=="cell_type_1")
two_cells <- df %>% dplyr::filter(cell_type=="cell_type_2")

res <- RANN::nn2(data = one_cells[,1:2], query = two_cells[,1:2], k = 1)
dist_to_nearest <- res$nn.dists[, 1]
two_cells$dist_nearest <- dist_to_nearest

res <- RANN::nn2(data = two_cells[,1:2], query = one_cells[,1:2], k = 1)
dist_to_nearest <- res$nn.dists[, 1]
one_cells$dist_nearest <- dist_to_nearest

coords_df = dplyr::left_join(df %>%
                               dplyr::mutate(score = cooccur_local_df$cooccur_local_cell_type_1_cell_type_2), 
                             rbind(one_cells,two_cells), by = c("x","y","cell_type"))
coef = cor.test(na.omit(coords_df)$score,na.omit(coords_df)$dist_nearest,method = "spearman",use = "pairwise.complete.obs")$estimate
pval = cor.test(na.omit(coords_df)$score,na.omit(coords_df)$dist_nearest,method = "spearman",use = "pairwise.complete.obs")$p.value

options(repr.plot.width=6, repr.plot.height=6)
ggplot(data = coords_df, aes(x = dist_nearest, y = score)) +
  geom_point_rast() +
  geom_smooth(
    method = "lm", se = FALSE, color = "blue"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Distance to nearest each other (µm)",
    y = paste("Co-occurence score"),
    title = paste0(test_type,"\nEst. interacted distance = ", distance_param), 
    subtitle = paste0(
      "Spearman's r = ", round(coef, 2),
      ", adj-p = ", signif(pval, 3)
    )
  )
```

``` output
Warning message in cor.test.default(na.omit(coords_df)$score, na.omit(coords_df)$dist_nearest, :
“Cannot compute exact p-value with ties”
Warning message in cor.test.default(na.omit(coords_df)$score, na.omit(coords_df)$dist_nearest, :
“Cannot compute exact p-value with ties”
`geom_smooth()` using formula = 'y ~ x'
“Removed 400 rows containing non-finite outside the scale range
(`stat_smooth()`).”
“Removed 400 rows containing missing values or values outside the scale range
(`geom_point()`).”
```

![](figures/SNA_tutorial_simulation/fig-06.png)

Check the ROC curve and AUC value to evaluate the performance of the
co-localization score

Check the distribution of the co-localization score by cell type

``` r

coords_df %>%
  ggplot(aes(x = cell_type, y = score, fill = cell_type)) +
  geom_violin(width = 1, alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.7) +
  #geom_jitter(width = 0.1, alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "Cell Type",
    y = "Co-occurence score",
    title = paste0(test_type,"\nEst. interacted distance = ", distance_param)
  )
```

![](figures/SNA_tutorial_simulation/fig-07.png)

``` r

sessionInfo()
```

``` output
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.6

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] doParallel_1.0.17     iterators_1.0.14      foreach_1.5.2        
 [4] ggrastr_1.0.2         ComplexHeatmap_2.18.0 circlize_0.4.15      
 [7] dplyr_1.1.4           magrittr_2.0.3        ggplot2_3.5.2        
[10] patchwork_1.2.0       spatialCooccur_0.1.0 

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3     jsonlite_1.8.8         shape_1.4.6           
  [4] magick_2.8.2           ggbeeswarm_0.7.2       spatstat.utils_3.1-2  
  [7] farver_2.1.1           GlobalOptions_0.1.2    vctrs_0.6.5           
 [10] ROCR_1.0-11            Cairo_1.6-2            spatstat.explore_3.2-6
 [13] base64enc_0.1-3        htmltools_0.5.8.1      sctransform_0.4.1     
 [16] parallelly_1.36.0      KernSmooth_2.23-22     htmlwidgets_1.6.4     
 [19] ica_1.0-3              plyr_1.8.9             plotly_4.10.4         
 [22] zoo_1.8-12             uuid_1.2-0             igraph_2.0.1.1        
 [25] mime_0.12              lifecycle_1.0.4        pkgconfig_2.0.3       
 [28] Matrix_1.6-5           R6_2.5.1               fastmap_1.1.1         
 [31] clue_0.3-65            fitdistrplus_1.1-11    future_1.33.1         
 [34] shiny_1.8.0            digest_0.6.35          colorspace_2.1-0      
 [37] S4Vectors_0.40.2       Seurat_5.0.1           tensor_1.5            
 [40] RSpectra_0.16-1        irlba_2.3.5.1          labeling_0.4.3        
 [43] progressr_0.14.0       fansi_1.0.6            spatstat.sparse_3.0-3 
 [46] mgcv_1.9-1             httr_1.4.7             polyclip_1.10-6       
 [49] abind_1.4-5            compiler_4.3.2         withr_3.0.0           
 [52] fastDummies_1.7.3      MASS_7.3-60.0.1        rjson_0.2.21          
 [55] tools_4.3.2            vipor_0.4.7            lmtest_0.9-40         
 [58] beeswarm_0.4.0         httpuv_1.6.14          future.apply_1.11.1   
 [61] goftest_1.2-3          glue_1.7.0             nlme_3.1-164          
 [64] promises_1.2.1         pbdZMQ_0.3-11          Rtsne_0.17            
 [67] cluster_2.1.6          reshape2_1.4.4         generics_0.1.3        
 [70] gtable_0.3.4           spatstat.data_3.0-4    tidyr_1.3.1           
 [73] data.table_1.15.0      sp_2.1-3               utf8_1.2.4            
 [76] BiocGenerics_0.48.1    spatstat.geom_3.2-8    RcppAnnoy_0.0.22      
 [79] ggrepel_0.9.5          RANN_2.6.1             pillar_1.9.0          
 [82] stringr_1.5.1          spam_2.10-0            IRdisplay_1.1         
 [85] RcppHNSW_0.6.0         later_1.3.2            splines_4.3.2         
 [88] moments_0.14.1         lattice_0.22-5         survival_3.5-7        
 [91] deldir_2.0-2           tidyselect_1.2.1       miniUI_0.1.1.1        
 [94] pbapply_1.7-2          gridExtra_2.3          IRanges_2.36.0        
 [97] scattermore_1.2        stats4_4.3.2           matrixStats_1.5.0     
[100] stringi_1.8.3          lazyeval_0.2.2         evaluate_0.23         
[103] codetools_0.2-19       tibble_3.2.1           cli_3.6.2             
[106] uwot_0.1.16            IRkernel_1.3.2         xtable_1.8-4          
[109] reticulate_1.39.0      repr_1.1.6             munsell_0.5.1         
[112] Rcpp_1.0.12            globals_0.16.2         spatstat.random_3.2-2 
[115] png_0.1-8              ellipsis_0.3.2         dotCall64_1.1-1       
[118] listenv_0.9.1          viridisLite_0.4.2      scales_1.3.0          
[121] ggridges_0.5.4         SeuratObject_5.0.2     leiden_0.4.3.1        
[124] purrr_1.0.2            crayon_1.5.2           GetoptLong_1.0.5      
[127] rlang_1.1.6            cowplot_1.1.3         
```
