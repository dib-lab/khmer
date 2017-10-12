# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=too-few-public-methods,no-init,missing-docstring
"""This is khmer; please see http://khmer.readthedocs.io/."""


from collections import namedtuple
from math import log
import json


from khmer._khmer import Read
from khmer._khmer import forward_hash
# tests/test_{functions,countgraph,counting_single}.py

from khmer._khmer import forward_hash_no_rc  # tests/test_functions.py

from khmer._khmer import reverse_hash  # tests/test_functions.py
# tests/counting_single.py

from khmer._khmer import hash_murmur3        # tests/test_functions.py
from khmer._khmer import hash_no_rc_murmur3  # tests/test_functions.py

from khmer._khmer import reverse_complement

from khmer._khmer import get_version_cpp as __version_cpp__
# tests/test_version.py

from khmer._khmer import ReadParser  # sandbox/to-casava-1.8-fastq.py
# tests/test_read_parsers.py,scripts/{filter-abund-single,load-graph}.py
# scripts/{abundance-dist-single,load-into-counting}.py

from khmer._khmer import FILETYPES

from khmer._oxli.graphs import (Counttable, QFCounttable, Nodetable,
                                CyclicCounttable,
                                SmallCounttable, Countgraph, SmallCountgraph,
                                Nodegraph)
from khmer._oxli.labeling import GraphLabels
from khmer._oxli.legacy_partitioning import SubsetPartition, PrePartitionInfo
from khmer._oxli.parsing import FastxParser
from khmer._oxli.readaligner import ReadAligner

from khmer._oxli.utils import get_n_primes_near_x, is_prime
import sys

from struct import pack, unpack

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


_buckets_per_byte = {
    # calculated by hand from settings in third-part/cqf/gqf.h
    'qfcounttable': 1 / 1.26,
    'countgraph': 1,
    'smallcountgraph': 2,
    'nodegraph': 8,
}


def extract_nodegraph_info(filename):
    """Open the given nodegraph file and return a tuple of information.

    Returns: the k-mer size, the table size, the number of tables, the version
    of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the nodegraph file to inspect
    """
    ksize = None
    n_tables = None
    table_size = None
    signature = None
    version = None
    ht_type = None
    occupied = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as nodegraph:
            signature, = unpack('4s', nodegraph.read(4))
            version, = unpack('B', nodegraph.read(1))
            ht_type, = unpack('B', nodegraph.read(1))
            ksize, = unpack('I', nodegraph.read(uint_size))
            n_tables, = unpack('B', nodegraph.read(uchar_size))
            occupied, = unpack('Q', nodegraph.read(ulonglong_size))
            table_size, = unpack('Q', nodegraph.read(ulonglong_size))
        if signature != b"OXLI":
            raise ValueError("Node graph '{}' is missing file type "
                             "signature".format(filename) + str(signature))
    except:
        raise ValueError("Node graph '{}' is corrupt ".format(filename))

    return ksize, round(table_size, -2), n_tables, version, ht_type, occupied


def extract_countgraph_info(filename):
    """Open the given countgraph file and return a tuple of information.

    Return: the k-mer size, the table size, the number of tables, the bigcount
    flag, the version of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the countgraph file to inspect
    """
    CgInfo = namedtuple("CgInfo", ['ksize', 'n_tables', 'table_size',
                                   'use_bigcount', 'version', 'ht_type',
                                   'n_occupied'])
    ksize = None
    n_tables = None
    table_size = None
    signature = None
    version = None
    ht_type = None
    use_bigcount = None
    occupied = None

    uint_size = len(pack('I', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as countgraph:
            signature, = unpack('4s', countgraph.read(4))
            version, = unpack('B', countgraph.read(1))
            ht_type, = unpack('B', countgraph.read(1))
            if ht_type != FILETYPES['SMALLCOUNT']:
                use_bigcount, = unpack('B', countgraph.read(1))
            else:
                use_bigcount = None
            ksize, = unpack('I', countgraph.read(uint_size))
            n_tables, = unpack('B', countgraph.read(1))
            occupied, = unpack('Q', countgraph.read(ulonglong_size))
            table_size, = unpack('Q', countgraph.read(ulonglong_size))
        if signature != b'OXLI':
            raise ValueError("Count graph file '{}' is missing file type "
                             "signature. ".format(filename) + str(signature))
    except:
        raise ValueError("Count graph file '{}' is corrupt ".format(filename))

    return CgInfo(ksize, n_tables, round(table_size, -2), use_bigcount,
                  version, ht_type, occupied)


def calc_expected_collisions(graph, force=False, max_false_pos=.2):
    """Do a quick & dirty expected collision rate calculation on a graph.

    Also check to see that collision rate is within threshold.

    Keyword argument:
    graph: the countgraph or nodegraph object to inspect
    """
    sizes = graph.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(graph.n_occupied())
    min_size = min(sizes)

    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht

    if fp_all > max_false_pos:
        print("**", file=sys.stderr)
        print("** ERROR: the graph structure is too small for ",
              file=sys.stderr)
        print("** this data set.  Increase data structure size",
              file=sys.stderr)
        print("** with --max_memory_usage/-M.", file=sys.stderr)
        print("**", file=sys.stderr)
        print("** Do not use these results!!", file=sys.stderr)
        print("**", file=sys.stderr)
        print("** (estimated false positive rate of %.3f;" % fp_all,
              file=sys.stderr, end=' ')
        print("max recommended %.3f)" % max_false_pos, file=sys.stderr)
        print("**", file=sys.stderr)

        if not force:
            sys.exit(1)

    return fp_all


from khmer._oxli.assembly import (LinearAssembler, SimpleLabeledAssembler,
                                  JunctionCountAssembler)
from khmer._oxli.hashset import HashSet
from khmer._oxli.hllcounter import HLLCounter
from khmer._oxli.labeling import GraphLabels
