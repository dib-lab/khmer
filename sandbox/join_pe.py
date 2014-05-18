#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
from screed.fasta import fasta_iter


def get_name(r):
    name = r['name'].split()[0]
    name = name.rsplit('/', 1)[0]
    return name


def _exhaust(h):
    while 1:
        yield h.next()


def join_pe(i1, i2, count_d):
    h1 = iter(i1)
    h2 = iter(i2)

    while 1:
        while 1:
            try:
                r1 = h1.next()
            except StopIteration:
                for r2 in _exhaust(h2):
                    yield r2
                return

            n1 = get_name(r1)
            if count_d[n1] != 1:
                break

            yield r1

        while 1:
            try:
                r2 = h2.next()
            except StopIteration:
                yield r1
                for r1 in _exhaust(h1):
                    yield r1
                return

            n2 = get_name(r2)
            if count_d[n2] != 1:
                break

            yield r2

        ###

        assert n1 == n2
        yield r1
        yield r2

filename1, filename2 = sys.argv[1:3]

name_count = {}
for n, record in enumerate(
        fasta_iter(open(filename1), parse_description=False)):
    if n % 10000 == 0:
        sys.stderr.write('...%d\n' % n)

    name = get_name(record)
    name_count[name] = 1

for n, record in enumerate(
        fasta_iter(open(filename2), parse_description=False)):
    if n % 10000 == 0:
        sys.stderr.write('...%d, x2\n' % n)

    name = get_name(record)
    name_count[name] = name_count.get(name, 0) + 1

i1 = fasta_iter(open(filename1), parse_description=False)
i2 = fasta_iter(open(filename2), parse_description=False)

for record in join_pe(i1, i2, name_count):
    print '>%s\n%s' % (record['name'], record['sequence'])

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
