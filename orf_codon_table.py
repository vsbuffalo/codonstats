from Bio import SeqIO
from Bio.Alphabet import IUPAC
import sys
import itertools
from collections import Counter

CODONS = [''.join(x) for x in (itertools.product(IUPAC.IUPACUnambiguousDNA().letters, repeat=3))]

# These are classification strings that add an additional column based
# on strings in header (usually about pseudogenes). Adjust to your
# headers.
TYPE_HEADER_STRINGS = dict(PSC="likely contains premature stop codon)",
                           MFPSC="likely majority frameshift and contains premature stop codon)",
                           MF="likely majority frameshift)")
NUCLEOTIDES = IUPAC.IUPACAmbiguousDNA().letters

def get_header_type(header):
    "Get the header type"
    for type_header, substring in TYPE_HEADER_STRINGS.items():
        if substring in header:
            return type_header
    return "NI" # no interference

def contains_masked(seq):
    """
    Return whether a SeqRecord contains soft masking.
    """
    lc_nuc = set(NUCLEOTIDES.lower())
    if str(seq).isupper():
        return "none"
    for letter in seq:
        if letter in lc_nuc:
            return "some"
    return "all"
    

def codon_frequency(seq):
    """
    Return dict of codon frequency.
    """
    empty = Counter(dict([(c, 0) for c in CODONS]))
    tmp = [(str(seq.upper()[pos:pos+3])) for
           pos in range(0, len(seq), 3)]
    tmp = filter(lambda x: len(x) == 3, tmp)
    return empty + Counter(tmp)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("usage: python orf_codon_table.py <orfs.fasta>\n"
                         "Assumes FASTA file is of ORFs\n\nerror: too few arguments.\n")
        sys.exit(1)
    with open(sys.argv[1]) as orf_file:
        sys.stdout.write('contig\tmasking\ttype\t' + '\t'.join(CODONS) + "\n")
        for record in SeqIO.parse(orf_file, "fasta"):
            header_type = get_header_type(record.description)
            cf = codon_frequency(record.seq)
            masked = contains_masked(record.seq)

            row = [record.id, masked, header_type]
            row.extend([str(cf[k]) for k in CODONS])
            sys.stdout.write("\t".join(row) + "\n")
