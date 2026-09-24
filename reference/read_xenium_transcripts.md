# Read a Xenium transcript table

\*\*Experimental.\*\* Read \`transcripts.parquet\` (or
\`transcripts.csv.gz\`) from a Xenium output folder into a plain
data.frame with columns \`x\`, \`y\`, \`z\`, \`gene\`, \`qv\` (plus
\`cell_id\` / \`overlaps_nucleus\` when present), ready for
\[bin_transcripts()\]. Handles older outputs in which \`feature_name\`
is stored as a binary column, filters on the Phred-scaled quality
\`qv\`, and drops negative-control and unassigned probes.

## Usage

``` r
read_xenium_transcripts(path, qv_min = 20, drop_controls = TRUE, bbox = NULL)
```

## Arguments

- path:

  A Xenium output folder or the path of the transcript file.

- qv_min:

  Minimum \`qv\` to keep (10x recommends 20). \`NULL\` keeps all.

- drop_controls:

  Drop features whose names start with \`NegControl\`, \`BLANK\`,
  \`Unassigned\`, \`Deprecated\`, \`antisense\` or \`Intergenic\`.

- bbox:

  Optional \`c(xmin, xmax, ymin, ymax)\` crop, applied while reading.

## Value

A data.frame.

## Examples

``` r
if (FALSE) { # \dontrun{
tx <- read_xenium_transcripts("Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs")
b <- bin_transcripts(tx, bin_size = 4)
} # }
```
