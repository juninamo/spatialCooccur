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

‘24 September, 2026’

## Spatial neighborhood analysis (SNA, cell type level analysis)

First generate dummy data where two cell types (cell_type_1 and
cell_type_2 in the “cell_type” column) are close to each other:
cell_type_2 forms a ring around a disk of cell_type_1 (concentric
circles), for a fraction `close_ratio` of the cells.

> **Note (v0.99.2).** Earlier versions shuffled the row and column
> labels independently in the permutation test. That inflated same-type
> z-scores to about +10 and biased different-type z-scores downwards,
> most strongly when there were few cell types (hence the old advice to
> use more than 10 cell types). The permutation now uses one shuffle for
> both, and the z-scores are calibrated under spatial randomness for any
> number of cell types (3 to 25 tested).

``` r

suppressPackageStartupMessages(suppressWarnings({
  if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(spatialCooccur)
  library(patchwork)
  library(ggplot2)
  library(magrittr)
  library(dplyr)
  library(circlize)
  library(ComplexHeatmap)
  library(ggrastr)
}))
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
| Clustercell_type_1 | 23.9462152 | 7.55240752 | -5.3958430 | -3.82373224 | -5.80635244 | -2.26689053 | -5.43668687 | -3.84599452 | -4.5108173 | -0.35758967 | -2.28902658 | -4.29912555 | -1.92571443 | -3.32436320 | -4.2660152 |
| Clustercell_type_2 | 11.6664396 | 10.54484214 | -3.0764086 | -0.82269048 | -3.23832605 | -3.49157735 | -3.07689551 | -4.17474364 | -2.6125001 | -3.52852118 | -0.83772671 | -3.27183402 | -1.61240564 | -2.55429496 | -3.8987714 |
| Clustercell_type_3 | -4.7903951 | -1.31431602 | -0.1065409 | 2.03924634 | 1.49391337 | 1.38742042 | 2.06612177 | 1.10913149 | 0.7275787 | 1.27762327 | 1.40551062 | 1.20882289 | -0.07352066 | 0.04911908 | 1.7309958 |
| Clustercell_type_4 | -2.5971542 | -0.01574833 | 0.7607068 | -0.45675781 | 0.27936360 | -0.39076603 | 0.46849941 | -0.19154561 | 0.1567979 | 0.38388104 | 1.32517014 | -0.66235798 | -1.11936001 | 0.74591614 | -0.7338555 |
| Clustercell_type_5 | -5.9018119 | -2.49389597 | 1.4220691 | 1.63744595 | 0.17649029 | 1.57788513 | 2.15393321 | 1.64587923 | 2.3365992 | 0.26511365 | -0.13822270 | 2.82392510 | 0.52710098 | 1.78246712 | 1.9532530 |
| Clustercell_type_6 | -0.7889704 | -2.09523368 | -0.2945250 | 0.39670465 | 0.26683303 | -0.37977428 | 0.28285391 | 0.39909189 | 1.5204003 | 1.01036955 | -1.30876962 | 1.17120012 | -0.73449688 | -0.39735893 | 2.7057134 |
| Clustercell_type_7 | -5.3189627 | -1.57800022 | 1.5167070 | 0.97343333 | 1.14744226 | 2.52227593 | 0.91148352 | 0.04390953 | 1.5291142 | 0.84523018 | -0.60576832 | 0.80226186 | -0.17905000 | 0.55670640 | 1.0836499 |
| Clustercell_type_8 | -4.2232122 | -3.47448344 | 0.8587543 | 0.92884826 | 1.50286681 | 0.78444010 | 0.38646304 | -0.16301631 | 0.8165288 | 0.80205192 | 1.02273357 | -0.51558702 | 1.44169053 | 2.46512023 | 1.6798632 |
| Clustercell_type_9 | -3.7111627 | -1.87235391 | -0.1871118 | 0.08593054 | 0.54578219 | 1.81785875 | 1.59159355 | 0.93051831 | 0.4194477 | 0.73241891 | -0.31794703 | 2.01237126 | 1.22422127 | -0.03826781 | 0.7692662 |
| Clustercell_type_10 | 1.2898659 | -2.95208987 | -0.6798872 | -0.17312173 | 0.04244105 | 1.65074606 | 0.96998256 | 0.39125354 | -0.6616165 | -1.55610914 | 0.04299035 | 0.48039933 | -1.07396824 | 0.44995656 | 1.3685586 |
| Clustercell_type_11 | -0.8340541 | 0.56242526 | 0.9748902 | 1.04756004 | -0.98921206 | -2.17068449 | -0.30223534 | 1.46187500 | -0.6870519 | 0.96698617 | -1.21993249 | 0.13879806 | 0.11199731 | 1.18853450 | -0.1425264 |
| Clustercell_type_12 | -3.4448219 | -1.66118506 | -0.2057989 | 0.05974900 | 0.28957481 | 0.08644692 | -0.03856619 | -0.65749684 | 0.9706269 | 0.39570032 | 0.32181942 | -0.02977462 | 0.80664099 | 0.53222064 | 0.7668497 |
| Clustercell_type_13 | -2.5814666 | -1.68473588 | 0.7105855 | 0.23894026 | 0.30891844 | -0.99346511 | 0.60728395 | 1.58949172 | 0.8372079 | -0.09399647 | 0.16496200 | 0.99344001 | 0.70798618 | 1.26178019 | 1.1170149 |
| Clustercell_type_14 | -2.2931761 | -0.76467865 | -0.3148586 | -0.17747055 | 1.31715287 | -1.45048262 | 0.18764913 | 1.10728429 | -1.3150995 | 0.08186561 | 0.06230033 | 0.63262575 | -0.68156380 | -0.62075647 | -0.4057012 |
| Clustercell_type_15 | -4.3679490 | -3.06176850 | -0.6248035 | 0.41598618 | 1.15980091 | 1.52799462 | 0.64019114 | 1.77480684 | 0.5590566 | 1.11072644 | 0.60358689 | 0.79449527 | 0.82226981 | 1.29512233 | -0.2911776 |

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

**z-score or log2 O/E?**
[`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
also returns `log2_oe`, log2(observed / expected). The z-score measures
the *evidence* for enrichment within this sample and grows with the
number of cells; `log2_oe` measures its *size*. Use the z-score within
one sample, as here, and `log2_oe` when comparing samples or groups (see
the case-control tutorial).

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
random_seeds <- sample(1000:9999, 20)
n_types_sim <- n_types

accuracy_df_all <- do.call(rbind, lapply(c(5, 10, 20, 30, 40, 50, 75, 100), function(distance_param) {
  do.call(rbind, lapply(random_seeds, function(seed_) {
    # a new tissue for every seed
    df <- generate_sim(close_ratio = close_ratio, n_types = n_types_sim, max_loc = max_loc,
                       n_cells = n_cells, test_type = "circle",
                       distance_param = distance_param, seed = seed_)
    z <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = neighbors.k_,
                          connectivity_key = "nn", transformation = TRUE,
                          n_perms = n_perm, seed = seed_, n_jobs = 1)$zscore
    dimnames(z) <- lapply(dimnames(z), function(v) gsub("^Cluster", "", v))
    data.frame(test_type = "circle", seed = seed_, distance_param = distance_param,
               zscore = z["cell_type_1", "cell_type_2"],        # planted pair
               zscore_false = z["cell_type_3", "cell_type_4"])  # unrelated pair
  }))
}))
accuracy_df_all$n_types <- n_types_sim
accuracy_df_all$neighbors.k_ <- neighbors.k_
head(accuracy_df_all)
```

|  | test_type | seed | distance_param | zscore | zscore_false | n_types | neighbors.k\_ |
|----|----|----|----|----|----|----|----|
|  | \<chr\> | \<int\> | \<dbl\> | \<dbl\> | \<dbl\> | \<dbl\> | \<dbl\> |
| 1 | circle | 8451 | 5 | 13.247010 | -1.89891730 | 15 | 30 |
| 2 | circle | 9015 | 5 | 13.753787 | 1.28154388 | 15 | 30 |
| 3 | circle | 8161 | 5 | 9.817993 | 0.60922574 | 15 | 30 |
| 4 | circle | 9085 | 5 | 14.703173 | 1.49370460 | 15 | 30 |
| 5 | circle | 8268 | 5 | 14.571439 | 0.09802873 | 15 | 30 |
| 6 | circle | 1622 | 5 | 16.173440 | -0.39799206 | 15 | 30 |

A data.frame: 6 × 7 {.table .dataframe}

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
```

``` output
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
 Min.   :0.001175                     
 1st Qu.:0.037670                     
 Median :0.103213                     
 Mean   :0.102000                     
 3rd Qu.:0.161757                     
 Max.   :0.231130                     
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
```

``` output
Warning message in cor.test.default(na.omit(coords_df)$score, na.omit(coords_df)$dist_nearest, :
“Cannot compute exact p-value with ties”
```

``` output
`geom_smooth()` using formula = 'y ~ x'
```

``` output
“Removed 400 rows containing non-finite values (`stat_smooth()`).”
```

``` output
“Removed 400 rows containing missing values (`geom_point()`).”
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
Running under: macOS 26.3.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] ggrastr_1.0.2         ComplexHeatmap_2.18.0 circlize_0.4.15      
[4] dplyr_1.1.4           magrittr_2.0.3        ggplot2_3.4.4        
[7] patchwork_1.1.3       spatialCooccur_0.99.2 testthat_3.2.1       

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.21       splines_4.3.2          later_1.3.2           
  [4] pbdZMQ_0.3-10          tibble_3.2.1           polyclip_1.10-6       
  [7] fastDummies_1.7.3      lifecycle_1.0.4        doParallel_1.0.17     
 [10] rprojroot_2.0.4        globals_0.16.2         lattice_0.21-9        
 [13] MASS_7.3-60            plotly_4.10.3          remotes_2.4.2.1       
 [16] httpuv_1.6.13          Seurat_5.2.1           sctransform_0.4.1     
 [19] spam_2.10-0            sp_2.1-2               sessioninfo_1.2.2     
 [22] pkgbuild_1.4.3         spatstat.sparse_3.1-0  reticulate_1.35.0     
 [25] cowplot_1.1.2          pbapply_1.7-2          RColorBrewer_1.1-3    
 [28] abind_1.4-5            pkgload_1.3.3          Rtsne_0.17            
 [31] purrr_1.0.2            BiocGenerics_0.48.1    IRanges_2.36.0        
 [34] S4Vectors_0.40.2       ggrepel_0.9.4          irlba_2.3.5.1         
 [37] listenv_0.9.0          spatstat.utils_3.1-2   moments_0.14.1        
 [40] goftest_1.2-3          RSpectra_0.16-1        spatstat.random_3.3-2 
 [43] fitdistrplus_1.1-11    parallelly_1.36.0      codetools_0.2-19      
 [46] tidyselect_1.2.0       shape_1.4.6            farver_2.1.1          
 [49] matrixStats_1.2.0      stats4_4.3.2           base64enc_0.1-3       
 [52] spatstat.explore_3.3-4 jsonlite_2.0.0         GetoptLong_1.0.5      
 [55] ellipsis_0.3.2         progressr_0.14.0       ggridges_0.5.5        
 [58] survival_3.5-7         iterators_1.0.14       foreach_1.5.2         
 [61] tools_4.3.2            ica_1.0-3              Rcpp_1.0.11           
 [64] glue_1.6.2             gridExtra_2.3          mgcv_1.9-0            
 [67] usethis_2.2.2          IRdisplay_1.1          withr_2.5.2           
 [70] fastmap_1.1.1          digest_0.6.33          R6_2.5.1              
 [73] mime_0.12              colorspace_2.1-0       scattermore_1.2       
 [76] Cairo_1.6-2            tensor_1.5             spatstat.data_3.1-4   
 [79] tidyr_1.3.0            generics_0.1.3         data.table_1.16.0     
 [82] httr_1.4.7             htmlwidgets_1.6.4      uwot_0.1.16           
 [85] pkgconfig_2.0.3        gtable_0.3.4           lmtest_0.9-40         
 [88] brio_1.1.4             htmltools_0.5.7        profvis_0.3.8         
 [91] dotCall64_1.1-1        clue_0.3-65            SeuratObject_5.0.2    
 [94] scales_1.3.0           png_0.1-8              spatstat.univar_3.1-2 
 [97] rstudioapi_0.15.0      reshape2_1.4.4         rjson_0.2.23          
[100] uuid_1.1-1             nlme_3.1-163           repr_1.1.6            
[103] cachem_1.0.8           zoo_1.8-12             GlobalOptions_0.1.2   
[106] stringr_1.5.1          KernSmooth_2.23-22     parallel_4.3.2        
[109] miniUI_0.1.1.1         vipor_0.4.7            desc_1.4.3            
[112] pillar_1.11.0          vctrs_0.6.5            RANN_2.6.1            
[115] urlchecker_1.0.1       promises_1.2.1         xtable_1.8-4          
[118] cluster_2.1.4          beeswarm_0.4.0         evaluate_0.23         
[121] magick_2.8.2           cli_3.6.2              compiler_4.3.2        
[124] rlang_1.1.2            crayon_1.5.2           future.apply_1.11.1   
[127] labeling_0.4.3         plyr_1.8.9             fs_1.6.3              
[130] ggbeeswarm_0.7.2       stringi_1.8.3          viridisLite_0.4.2     
[133] deldir_2.0-2           munsell_0.5.0          lazyeval_0.2.2        
[136] devtools_2.4.5         spatstat.geom_3.3-5    Matrix_1.6-5          
[139] IRkernel_1.3.2         RcppHNSW_0.5.0         future_1.33.1         
[142] shiny_1.8.0            ROCR_1.0-11            igraph_1.6.0          
[145] memoise_2.0.1         
```
