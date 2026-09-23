# Algorithm reference: core spatialCooccur functions

This vignette documents the mathematics behind each core function in
**spatialCooccur**. The formulas mirror the implementation in
`R/spatialCooccur_functions.R` line-by-line, not a stylized version.

``` r

library(spatialCooccur)
#> Loading spatialCooccur v0.99.2: An R package for analyzing spatial co-occurrence.
#> To cite this package in publications, please use:
#>   Inamo J, et al. (2026). Spatial transcriptomics reveals  immune-stromal crosstalk within the synovium of patients with  juvenile idiopathic arthritis. JCI Insight. 11(1):e198074.  doi:10.1172/jci.insight.198074
#> Developed by: Jun Inamo <juninamo@keio.jp>
```

## Notation

Throughout, a sample is a finite set of cells

``` math
\mathcal{V} = \{v_1, \dots, v_N\}, \quad
\mathbf{x}_v \in \mathbb{R}^{2}, \quad
c(v) \in \mathcal{C} = \{c_1, \dots, c_K\},
```

where $`\mathbf{x}_v`$ is the 2-D centroid and $`c(v)`$ the discrete
cluster label (e.g. cell type). For each function we additionally use:

- $`A \in \{0, 1\}^{N \times N}`$ — adjacency of a
  $`k`$-nearest-neighbor (kNN) graph on $`\{\mathbf{x}_v\}`$, with
  $`A_{uv} = 1`$ iff $`v`$ is one of the $`k`$ nearest neighbors of
  $`u`$. Built with `Seurat::FindNeighbors(..., k.param = k)`.
- $`d_v = \sum_u A_{uv}`$ — the (in-)degree of $`v`$ in the kNN graph.
- $`N_r(v) = \{u \in \mathcal{V}\setminus\{v\} : \lVert \mathbf{x}_u -
  \mathbf{x}_v \rVert \le r\}`$ — radius-$`r`$ neighborhood, computed
  with `RANN::nn2(..., searchtype = "radius", radius = r, k = k_\max)`.
- $`\mathcal{V}_c = \{v : c(v) = c\}`$ — the set of cells in cluster
  $`c`$.

## `generate_sim()` — synthetic spatial layouts

[`generate_sim()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim.md)
produces a single sample under one of three generative schemes. The user
supplies `close_ratio` $`\rho`$, `distance_param` $`\Delta`$, `n_cells`
$`N`$, `n_types` $`K`$, and `max_loc` $`L`$ (the side length of the
bounding box).

In every case the *non-special* cells are drawn uniformly from
$`\mathcal{U}\bigl([0, L]^2\bigr)`$ and the cluster labels are sampled
i.i.d. from a uniform distribution over $`\{c_1, \dots, c_K\}`$. What
differs between `test_type`s is how `cell_type_1` and `cell_type_2` are
positioned.

### `test_type = "distribute"`

Start from $`N`$ i.i.d. samples in $`[0, L]^2`$ with i.i.d. cluster
labels. Let $`I_k = \{v : c(v) = c_k\}`$ denote the index set of cluster
$`k`$. A fraction $`\rho`$ of `cell_type_2` cells, of size

``` math
n_\text{pair} = \mathrm{round}\bigl(\rho \cdot |I_2|\bigr),
```

is randomly chosen and *moved* to be near a uniformly random
`cell_type_1` cell:

``` math
\theta_m \stackrel{\text{iid}}{\sim} \mathcal{U}(0, 2\pi),
\quad
\boldsymbol{\eta}_m \stackrel{\text{iid}}{\sim} \mathcal{N}\!\left(\mathbf{0},\, (\Delta/5)^2 I_2\right),
```

``` math
\mathbf{x}_{v_m^{(2)}} \;\leftarrow\;
  \mathbf{x}_{v_m^{(1)}} \;+\; \Delta \binom{\cos \theta_m}{\sin \theta_m} \;+\; \boldsymbol{\eta}_m,
\quad m = 1, \dots, n_\text{pair},
```

where $`(v_m^{(1)}, v_m^{(2)})`$ are randomly paired indices from
$`I_1`$ and $`I_2`$. So $`\rho`$ controls *how many* type-2 cells are
seeded next to a type-1 cell, and $`\Delta`$ controls *how close* the
pair is.

### `test_type = "circle"`

[`generate_sim()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim.md)
sets the disk centre at $`\mathbf{c} = (L/2, L/2)`$ and the radius at
$`R = \lceil L/6 \rceil`$. Of $`n_1 = \lceil N/K \rceil`$`cell_type_1`
cells, a fraction $`\rho`$ is placed inside / on the disk:

``` math
\begin{aligned}
n_1^{\text{close}} &= \mathrm{round}(\rho \, n_1), \\
\text{outer ring:} \quad
\mathbf{x}_v &= \mathbf{c} + R\,(\cos\theta, \sin\theta),
\quad \theta \in \text{uniform grid on } [0, 2\pi], \\
\text{inner disk:} \quad
\mathbf{x}_v &= \mathbf{c} + \rho_v (\cos\theta_v, \sin\theta_v),
\quad \rho_v = \sqrt{U_v}\,R, \;
U_v \sim \mathcal{U}(0, 1), \;
\theta_v \sim \mathcal{U}(0, 2\pi),
\end{aligned}
```

with the remaining $`n_1 - n_1^\text{close}`$ cells dropped uniformly in
$`[0, L]^2`$. Of $`n_2 = \lceil N/K \rceil`$`cell_type_2` cells,
$`n_2^\text{close} = \mathrm{round}(\rho \, n_2)`$ are placed on the
outer ring $`R + \Delta`$:

``` math
\mathbf{x}_v = \mathbf{c} + (R + \Delta)\,(\cos\theta, \sin\theta) + \boldsymbol{\eta},
\quad \boldsymbol{\eta} \sim \mathcal{N}(\mathbf{0}, 5^2 I_2).
```

The remainder of `cell_type_2` is uniform, as are the other $`K - 2`$
clusters. The result is the familiar “type-1 disk + type-2 ring” layout:
the two clusters co-occur at the scale of the whole disk, but a
small-$`k`$ local neighborhood graph mostly sees same-cluster neighbors.

### `test_type = "line"`

Two thin horizontal layers separated vertically by $`\Delta`$. Layer 0
contains `cell_type_1` along an evenly spaced $`x`$-axis grid with
Gaussian jitter; layer 1 contains `cell_type_2` at $`y`$-position
shifted by $`\Delta`$. A fraction $`\rho`$ of layer-1 cells are
*snapped* to the $`(x, y + \Delta)`$ position of a randomly chosen
layer-0 cell, so that $`\rho`$ again controls how many `cell_type_2`
cells are tightly paired with a `cell_type_1` cell.

## `nhood_enrichment()` — neighborhood enrichment z-score

This is the heaviest of the core functions and the most natural to
formalize.

#### Step 1 — build the kNN adjacency

Construct the kNN graph with
`Seurat::FindNeighbors(coords, k.param = k)`. Let
$`A \in \{0, 1\}^{N \times N}`$ be the resulting adjacency
(`neighbors$nn`) or the shared-NN matrix (`neighbors$snn`), depending on
`connectivity_key`.

#### Step 2 — degree normalize (optional)

When `transformation = TRUE`, divide each column by the (column) degree
plus one to keep highly-connected cells from dominating:

``` math
\tilde A_{uv} \;=\; \frac{A_{uv}}{d_v + 1},
\qquad
d_v = \sum_{u} A_{uv}.
```

Otherwise $`\tilde A = A`$ and the function counts edges as binary.

#### Step 3 — cluster-pair count matrix

For every ordered cluster pair $`(c_i, c_j)`$, sum entries of
$`\tilde A`$ whose row cell belongs to $`c_i`$ and column cell to
$`c_j`$:

``` math
C_{ij} \;=\; \sum_{u \in \mathcal{V}_{c_i}} \sum_{v \in \mathcal{V}_{c_j}} \tilde A_{uv}.
```

(With `transformation = FALSE` the inner term is $`\mathbb{1}\{A_{uv} =
1\}`$ instead — i.e. raw edge counts.)

This is the observed value. Implemented in
[`compute_count()`](https://juninamo.github.io/spatialCooccur/reference/compute_count.md).

#### Step 4 — permutation null

[`permute_clusters()`](https://juninamo.github.io/spatialCooccur/reference/permute_clusters.md)
shuffles the cluster vector $`\mathbf{c} = (c(v_1),
\dots, c(v_N))`$ twice — *independently* for rows and columns — and
recomputes the count:

``` math
C^{(b)}_{ij} \;=\; \sum_{u \in \mathcal{V}^{(\pi^{(b)}_r)}_{c_i}}
                   \sum_{v \in \mathcal{V}^{(\pi^{(b)}_c)}_{c_j}}
                   \tilde A_{uv},
\quad b = 1, \dots, B,
```

where $`\pi^{(b)}_r, \pi^{(b)}_c`$ are independent uniform random
permutations of $`\{1, \dots, N\}`$ and $`\mathcal{V}^{(\pi)}_c`$ is the
set of cells whose *permuted* label equals $`c`$. The permutations
destroy any association between cluster identity and spatial position,
giving the null distribution of $`C_{ij}`$.

#### Step 5 — z-score

Let $`\hat\mu_{ij}, \hat\sigma_{ij}`$ be the empirical mean and standard
deviation of $`\{C^{(b)}_{ij}\}_{b=1}^B`$. The reported enrichment
z-score is

``` math
\boxed{\;
z_{ij} \;=\; \frac{C_{ij} - \hat\mu_{ij}}{\hat\sigma_{ij}}\;.
}
```

Positive $`z_{ij}`$ means cluster $`c_i`$ and cluster $`c_j`$ co-occur
more than expected under shuffling; negative means segregation.

In code:

``` r

res <- nhood_enrichment(
  df,
  cluster_key  = "cell_type",
  neighbors.k  = 30,
  connectivity_key = "nn",
  transformation = TRUE,
  n_perms = 100
)
res$count    # the C_{ij} matrix
res$zscore   # the z_{ij} matrix
```

## `calc_co_occurrence_for_radius()` and `compute_co_occurrence_ratio()`

Where
[`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)
works on a fixed-$`k`$ graph, these two functions work on a
fixed-*radius* neighborhood and report a conditional-probability ratio
instead of a z-score.

#### Step 1 — radius-based neighbor counts

For each cell $`v`$, find its radius neighborhood $`N_r(v)`$ via
[`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)
(with cap $`k`$). Tally pairs by cluster:

``` math
C_{ij} \;=\; \sum_{v \in \mathcal{V}_{c_i}}\; \bigl|N_r(v) \cap \mathcal{V}_{c_j}\bigr|.
```

Note that the same physical pair $`(v, u)`$ contributes once to
$`C_{c(v), c(u)}`$ and once to $`C_{c(u), c(v)}`$, so $`C`$ is symmetric
when accumulated over all cells. The matrix is stored under
`co_occur_count`.

#### Step 2 — enrichment ratio

[`compute_co_occurrence_ratio()`](https://juninamo.github.io/spatialCooccur/reference/compute_co_occurrence_ratio.md)
then converts $`C`$ to a ratio of empirical probabilities (Lifshitz /
Giotto convention):

``` math
P(j \mid i) \;=\; \frac{C_{ij}}{\sum_k C_{ik}},
\qquad
P(j) \;=\; \frac{\sum_k C_{kj}}{\sum_{k, l} C_{kl}},
```

``` math
\boxed{\;
r_{ij} \;=\; \frac{P(j \mid i)}{P(j)}\;.
}
```

$`r_{ij} > 1`$ means a cell in cluster $`c_i`$ is more likely to have a
cluster-$`c_j`$ neighbor than chance predicts from the overall pair
composition. Unlike the z-score in
§[`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md),
this ratio is deterministic — no permutation is run.

## `cooccur_local()` — per-cell local co-occurrence

Where the previous two functions produce a cluster-pair matrix, this one
produces a *per-cell score* for a single pair $`(c_x, c_y)`$ chosen by
the caller. It is meant to highlight individual cells that sit at the
interface of two cluster populations.

#### Step 1 — local presence indicator

Take the radius-$`r`$ neighborhood $`N_r(v)`$ (with neighbor cap $`k`$).
The raw local co-occurrence score is

``` math
s_v^{(0)} \;=\; \mathbb{1}\!\left\{
  N_r(v) \cap \mathcal{V}_{c_x} \ne \emptyset
  \;\wedge\;
  N_r(v) \cap \mathcal{V}_{c_y} \ne \emptyset
\right\}.
```

So $`s_v^{(0)} = 1`$ exactly when at least one cell of each target
cluster appears within distance $`r`$ of $`v`$ (and $`v`$ itself is
excluded). The score is symmetric in $`(c_x, c_y)`$.

#### Step 2 — single graph-diffusion step

The raw indicator is sparse.
[`cooccur_local()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.md)
smooths it once over the kNN graph $`A`$ (same construction as in
[`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md)).
With $`D = \mathrm{diag}(d_v + 1)`$:

``` math
\boxed{\;
s_v^{(1)} \;=\; \frac{s_v^{(0)}}{d_v + 1} \;+\; \sum_{u} A_{vu} \, \frac{s_u^{(0)}}{d_u + 1},
\quad\text{i.e.}\quad
\mathbf{s}^{(1)} \;=\; (A + I)\, D^{-1}\, \mathbf{s}^{(0)}.
\;}
```

Concretely, every cell shares an equal portion of its own indicator with
each of its kNN neighbors plus itself. The output column is named
`cooccur_local_<cluster_x>_<cluster_y>`.

With `maxnsteps > 1` the step is iterated, each step starting from the
previous result:
$`\mathbf{s}^{(t)} = (A + I) D^{-1} \mathbf{s}^{(t-1)}`$. The total
$`\sum_v s_v`$ is conserved. Iteration stops early when the kurtosis of
the scores decreases by less than 3 between consecutive steps (checked
after step 3). Because $`\mathbf{s}^{(0)}`$ is a 0/1 indicator, its
kurtosis usually *increases* under smoothing, so in practice diffusion
stops after 4 steps whenever `maxnsteps >= 4`. `maxnsteps = 0` returns
the raw indicator. (Versions \<= 0.99.1 recomputed every step from
$`\mathbf{s}^{(0)}`$, so any `maxnsteps >= 1` gave the single step
above.)

## `search_interaction_spot()` — connected-component spots

This function turns a radius-neighborhood graph into discrete
“interaction spots” — connected groups of cells that are likely to be
acting on each other in a tissue.

#### Step 1 — radius-neighbor graph

Restrict to a user-provided cell ID list $`\mathcal{V}^\star \subseteq
\mathcal{V}`$. For each $`v \in \mathcal{V}^\star`$ compute $`N_r(v)`$
via [`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)
(radius search, neighbor cap $`k`$). Build an undirected graph
$`G^\star = (\mathcal{V}^\star, E^\star)`$ with

``` math
E^\star \;=\; \bigl\{\{v, u\} : v \in \mathcal{V}^\star,\; u \in N_r(v) \cap \mathcal{V}^\star\bigr\}.
```

#### Step 2 — weak connected components

Run `igraph::components(G^\star, mode = "weak")`. Each component
$`\Gamma_q \subseteq \mathcal{V}^\star`$ is a candidate spot, and every
cell $`v`$ receives an integer label $`q(v)`$.

#### Step 3 — keep target-cluster, size-$`n_\min`$ spots

Define the bounding box of component $`q`$ from its cells’ coordinates,

``` math
\bigl[x_q^{\min}, x_q^{\max}\bigr] \times \bigl[y_q^{\min}, y_q^{\max}\bigr],
\quad
n_q = |\Gamma_q|,
```

and keep cell $`v`$ in the output iff

``` math
c(v) \in \mathcal{T} \quad\text{and}\quad n_{q(v)} > n_\min,
```

where $`\mathcal{T}`$ is the `target_cluster` argument. The function
returns a tidy data.frame of surviving cells, one row per cell, with
columns `cluster_id`, `x_min`, `x_max`, `y_min`, `y_max`, and
`n_all_cells`.

## Summary

| Function | Per-cell or pairwise | Local geometry | Statistic | Test |
|----|----|----|----|----|
| [`nhood_enrichment()`](https://juninamo.github.io/spatialCooccur/reference/nhood_enrichment.md) | Cluster pair $`(i, j)`$ | $`k`$-NN graph $`A`$ | $`z`$-score of $`C_{ij}`$ | Permutation null |
| [`calc_co_occurrence_for_radius()`](https://juninamo.github.io/spatialCooccur/reference/calc_co_occurrence_for_radius.md) / [`compute_co_occurrence_ratio()`](https://juninamo.github.io/spatialCooccur/reference/compute_co_occurrence_ratio.md) | Cluster pair $`(i, j)`$ | Radius-$`r`$ neighbors | Ratio $`P(j\mid i)/P(j)`$ | Deterministic |
| [`cooccur_local()`](https://juninamo.github.io/spatialCooccur/reference/cooccur_local.md) | Per cell, for fixed $`(c_x, c_y)`$ | Radius-$`r`$ + kNN diffusion | Real-valued score $`s_v^{(1)}`$ | — |
| [`search_interaction_spot()`](https://juninamo.github.io/spatialCooccur/reference/search_interaction_spot.md) | Per cell, group label | Radius-$`r`$ + connected comp. | Component size | — |
| [`generate_sim()`](https://juninamo.github.io/spatialCooccur/reference/generate_sim.md) | — | Generative | Spatial layout | — |

The disease-group extension (`*_per_sample()`,
[`compare_groups()`](https://juninamo.github.io/spatialCooccur/reference/compare_groups.md))
sits on top of these by running them per sample and testing the
resulting per-sample scalars between groups. See the [disease comparison
vignette](https://juninamo.github.io/spatialCooccur/articles/disease_comparison.md)
for the group-comparison side of the story.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] spatialCooccur_0.99.2
#> 
#> loaded via a namespace (and not attached):
#>   [1] deldir_2.0-4           pbapply_1.7-5          gridExtra_2.3.1       
#>   [4] rlang_1.3.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
#>   [7] otel_0.2.0             spatstat.geom_3.8-3    matrixStats_1.5.0     
#>  [10] ggridges_0.5.7         compiler_4.6.1         png_0.1-9             
#>  [13] systemfonts_1.3.2      vctrs_0.7.3            reshape2_1.4.5        
#>  [16] stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
#>  [19] promises_1.5.0         rmarkdown_2.32         ragg_1.5.2            
#>  [22] purrr_1.2.2            xfun_0.61              cachem_1.1.0          
#>  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8           
#>  [28] spatstat.utils_3.2-5   irlba_2.3.7            parallel_4.6.1        
#>  [31] cluster_2.1.8.2        R6_2.6.1               ica_1.0-3             
#>  [34] spatstat.data_3.1-9    bslib_0.12.0           stringi_1.8.9         
#>  [37] RColorBrewer_1.1-3     reticulate_1.47.0      spatstat.univar_3.2-0 
#>  [40] parallelly_1.48.0      lmtest_0.9-40          jquerylib_0.1.4       
#>  [43] scattermore_1.2        Rcpp_1.1.2             knitr_1.52            
#>  [46] tensor_1.5.1           future.apply_1.20.2    zoo_1.9-0             
#>  [49] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-5          
#>  [52] splines_4.6.1          igraph_2.3.3           tidyselect_1.2.1      
#>  [55] abind_1.4-8            yaml_2.3.12            spatstat.random_3.5-2 
#>  [58] spatstat.explore_3.8-3 codetools_0.2-20       miniUI_0.1.2          
#>  [61] listenv_1.0.0          lattice_0.22-9         tibble_3.3.1          
#>  [64] plyr_1.8.9             shiny_1.14.0           S7_0.2.2              
#>  [67] ROCR_1.0-12            evaluate_1.0.5         Rtsne_0.17            
#>  [70] future_1.75.0          fastDummies_1.7.6      desc_1.4.3            
#>  [73] survival_3.8-6         polyclip_1.10-7        fitdistrplus_1.2-6    
#>  [76] pillar_1.11.1          Seurat_5.5.1           KernSmooth_2.23-26    
#>  [79] plotly_4.12.1          generics_0.1.4         RcppHNSW_0.7.0        
#>  [82] sp_2.2-3               ggplot2_4.0.3          scales_1.4.0          
#>  [85] globals_0.19.1         xtable_1.8-8           glue_1.8.1            
#>  [88] tools_4.6.1            data.table_1.18.6.1    RSpectra_0.16-2       
#>  [91] RANN_2.6.3             fs_2.1.0               dotCall64_1.2         
#>  [94] cowplot_1.2.0          grid_4.6.1             tidyr_1.3.2           
#>  [97] nlme_3.1-169           patchwork_1.3.2        cli_3.6.6             
#> [100] spatstat.sparse_3.2-0  textshaping_1.0.5      spam_2.11-4           
#> [103] viridisLite_0.4.3      dplyr_1.2.1            uwot_0.2.5            
#> [106] gtable_0.3.6           sass_0.4.10            digest_0.6.39         
#> [109] progressr_1.0.0        ggrepel_0.9.8          htmlwidgets_1.6.4     
#> [112] SeuratObject_5.4.0     farver_2.1.2           htmltools_0.5.9       
#> [115] pkgdown_2.2.1          lifecycle_1.0.5        httr_1.4.9            
#> [118] mime_0.13              MASS_7.3-65
```
