#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
# Convenience functions for performing common argument-checking tasks in
# scripts.


def print_error(msg):
    """
        Prints the given message to 'stderr'.
    """

    import sys

    print >>sys.stderr, msg

def iter_kmers(seq, k):
    n = (len(seq) - k) + 1
    for i in xrange(n):
        yield seq[i:i+k]

def comp(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    else:
        return 'G'
   
def reverse_comp(seq):
    revc = [comp(b) for b in seq]
    for i in xrange(len(revc)/2):
        j = len(seq)-i-1
        if i < j:
            tmp = revc[i]
            revc[i] = revc[j]
            revc[j] = tmp
    return ''.join(revc)

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
