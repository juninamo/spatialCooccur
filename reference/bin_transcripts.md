# Bin transcript coordinates on a square grid

\*\*Experimental.\*\* Count transcripts per gene in square bins, without
cell segmentation. Every grid bin in the bounding box is kept (so that
the grid can be used with FFTs), and an \`in_tissue\` flag marks bins
close to at least one transcript.

## Usage

``` r
bin_transcripts(
  tx,
  bin_size = 4,
  x_col = "x",
  y_col = "y",
  gene_col = "gene",
  genes = NULL,
  qv_col = NULL,
  qv_min = 20,
  bbox = NULL,
  tissue_radius = 2
)
```

## Arguments

- tx:

  A data.frame of transcripts, e.g. from Xenium \`transcripts.parquet\`
  or \[simulate_transcripts()\].

- bin_size:

  Bin side length (same unit as the coordinates, usually micrometres).

- x_col, y_col, gene_col:

  Column names of coordinates and gene.

- genes:

  Genes to keep (default: all). Control probes such as \`NegControl\*\`
  / \`BLANK\*\` should be dropped here.

- qv_col, qv_min:

  Optional quality column and minimum value (Xenium: \`qv\`, 20).

- bbox:

  Optional bounding box \`c(xmin, xmax, ymin, ymax)\`.

- tissue_radius:

  Bins within this many bins of a transcript are flagged \`in_tissue\`.
  \`Inf\` flags every bin.

## Value

An object of class \`binned_transcripts\`: a list with \`counts\`
(sparse bins x genes matrix), \`coords\` (bin centres and
\`in_tissue\`), \`grid\` (\`nx\`, \`ny\`, \`bin_size\`, \`xmin\`,
\`ymin\`) and \`genes\`.

## Examples

``` r
tx <- simulate_transcripts(size = 100, rate = 0.01)
b <- bin_transcripts(tx, bin_size = 5)
dim(b$counts)
#> [1] 400  15
```
