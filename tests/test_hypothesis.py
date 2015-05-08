from hypothesis import given, assume
from hypothesis.specifiers import strings
from nose.plugins.attrib import attr

from khmer import reverse_hash, forward_hash, forward_hash_no_rc

#from hypothesis import Settings, Verbosity
#Settings.default.verbosity = Verbosity.verbose
TRANSLATE = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}


@attr('hypothesis')
@given(strings("ACGT"))
def test_forward_hash(kmer):
    ksize = len(kmer)
    assume(ksize > 0)
    assume(ksize < 32)

    rh = reverse_hash(forward_hash(kmer, ksize), ksize)
    assert rh == kmer or rh == "".join(TRANSLATE[c] for c in kmer[::-1])


@attr('hypothesis')
@given(strings("ACGT"))
def test_forward_hash_no_rc(kmer):
    ksize = len(kmer)
    assume(ksize > 0)
    assume(ksize < 32)

    rh = reverse_hash(forward_hash_no_rc(kmer, ksize), ksize)
    assert rh == kmer
