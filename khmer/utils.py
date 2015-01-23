#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
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

def write_record(record, fp):
    """
        Writes output sequence and returns str.
    """
    if hasattr(record, 'accuracy'):
        fp.write(
            '@{name}\n{seq}\n'
            '+\n{acc}\n'.format(name=record.name,
                               seq=record.sequence,
                               acc=record.accuracy))
    else:
       fp.write(
            '>{name}\n{seq}\n'.format(name=record.name,
                                      seq=record.sequence))

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
