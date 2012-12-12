from nose.tools import assert_equal
from orf_codon_table import codon_frequency
from collections import Counter
from operator import itemgetter

def test_codon_frequency():
    seq1 = "ATGGGACGACCCCCC"
    seq2 = "ATGGGACGACCCCCCT" # should not change counts
    seq1_counter = sorted(Counter(ATG=1, GGA=1, CGA=1, CCC=2).items(), key=itemgetter(0))
    assert_equal(sorted(codon_frequency(seq1).items(), key=itemgetter(0)), seq1_counter)
    assert_equal(sorted(codon_frequency(seq2).items(), key=itemgetter(0)), seq1_counter)



