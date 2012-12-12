# Codon Analysis

These are some simple programs for doing codon frequency analysis
using MDS plots. Most of this is currently geared around trying to do
this for results from `findorf`, so it may belong in `findorf/contrib`
more than a seperate module.

## Operation

    python orf_codon_table.py orfs.fasta > orf_table.tsv

    Rscript codon_stats.R orf_table.tsv <sample-size>

Where the `<sample-size>` corresponds to the number of ORFs to sample
from the psuedogene (majority frameshift or internal stop codon)
set. This measure uses equal numbers from both groups, even though
there are far more non-pseudogenes. This equal sampling approach is
just for graphical clarity. This assumes that there are more
pseudogenes than non-pseudogenes.

## Todo

 - Currently the "no interference" group is anything without a
   majority frameshift or internal stop codon (or both). We may wish
   to subset by full length ORFs only, in case partial ORFs are more
   likely to be frameshift or internal stop codon containing, and are
   not classified as such.