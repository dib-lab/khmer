#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
# Convenience functions for performing common argument-checking tasks in
# scripts.


def print_error(msg):
    """
        Prints the given message to 'stderr'.
    """

    import sys

    print >>sys.stderr, msg


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
