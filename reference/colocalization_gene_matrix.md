# Gene-by-gene co-localization of transcripts

\*\*Experimental.\*\* For genes \\a, b\\, counts transcript pairs within
distance \`radius\`, \$\$P\_{ab} = \sum_u Y_a(u)\\(K_r \* Y_b)(u),\$\$
where \\Y\\ are binned counts and \\K_r\\ is a disc of radius \`radius\`
(FFT convolution, one per gene), and compares it with its expectation
when gene labels are shuffled over the fixed transcript positions,
\$\$E\_{ab} = \frac{n_a n_b}{N(N-1)} P\_{\mathrm{all}},\$\$ with \\n\\
the transcript counts, \\N\\ their total and \\P\_{\mathrm{all}}\\ the
number of all transcript pairs within the radius (self pairs excluded).
\`log2((P + c) / (E + c))\` is 0 without co-localization; cellularity
and tissue shape cancel, as in the relative pair correlation.

## Usage

``` r
colocalization_gene_matrix(
  binned,
  genes = NULL,
  radius = 20,
  min_count = 50,
  top_n = NULL,
  pseudocount = 1
)
```

## Arguments

- binned:

  Output of \[bin_transcripts()\].

- genes:

  Genes to use (default: all genes with at least \`min_count\`
  transcripts).

- radius:

  Distance in um within which transcripts count as a pair.

- min_count:

  Minimum number of transcripts per gene.

- top_n:

  Keep at most this many genes (the most abundant).

- pseudocount:

  Added to observed and expected pair counts.

## Value

A symmetric gene x gene matrix of log2 O/E with attributes
\`n_transcripts\` (per gene) and \`radius\`.

## Examples

``` r
tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
M <- colocalization_gene_matrix(b, radius = 12)
round(M[1:4, 1:4], 2)
#>      A_1  A_2  A_3  A_4
#> A_1 0.80 0.82 0.78 0.82
#> A_2 0.82 0.86 0.81 0.84
#> A_3 0.78 0.81 0.76 0.80
#> A_4 0.82 0.84 0.80 0.83
```
