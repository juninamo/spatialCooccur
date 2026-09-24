# Cluster genes into co-localization modules

\*\*Experimental.\*\* Hierarchical clustering of the gene x gene log2
O/E matrix from \[colocalization_gene_matrix()\]: genes whose
transcripts lie near one another (high mutual O/E) end up in the same
module. The distance between genes is \`max(M) - M\` with average
linkage.

## Usage

``` r
colocalization_modules(M, n_modules = NULL, min_oe = 0.5, min_size = 3)
```

## Arguments

- M:

  Output of \[colocalization_gene_matrix()\].

- n_modules:

  Number of modules (\`cutree(k = )\`). If \`NULL\`, the tree is cut at
  height \`max(M) - min_oe\`, i.e. genes join a module when their
  average mutual log2 O/E exceeds \`min_oe\`.

- min_oe:

  Threshold used when \`n_modules\` is \`NULL\`.

- min_size:

  Modules with fewer genes are labelled \`NA\`.

## Value

A list with \`modules\` (data.frame: gene, module, connectivity = mean
log2 O/E with the other genes of its module), \`summary\` (data.frame:
module, size, mean_oe, top_genes) and the \`tree\`.

## Examples

``` r
tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
mods <- colocalization_modules(colocalization_gene_matrix(b, radius = 12), n_modules = 3)
mods$summary
#>   module size     mean_oe               top_genes
#> 1     M1    5 0.785425551 A_2, A_4, A_1, A_3, A_5
#> 2     M2    5 0.324309289 B_2, B_3, B_5, B_4, B_1
#> 3     M3    5 0.001650131 C_5, C_1, C_4, C_2, C_3
```
