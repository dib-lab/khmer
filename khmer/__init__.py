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


from math import log
import json


from khmer._khmer import Read
# tests/test_version.py

from khmer._khmer import ReadParser  # sandbox/to-casava-1.8-fastq.py
# tests/test_read_parsers.py,scripts/{filter-abund-single,load-graph}.py
# scripts/{abundance-dist-single,load-into-counting}.py

from khmer._oxli.assembly import (LinearAssembler, SimpleLabeledAssembler,
                                  JunctionCountAssembler)

from khmer._oxli.graphs import (Counttable, QFCounttable, Nodetable,
                                SmallCounttable, Countgraph, SmallCountgraph,
                                Nodegraph, _buckets_per_byte)

from khmer._oxli.hashing import (forward_hash, forward_hash_no_rc,
                                 reverse_hash, hash_murmur3,
                                 hash_no_rc_murmur3,
                                 reverse_complement)

from khmer._oxli.hashset import HashSet

from khmer._oxli.hllcounter import HLLCounter

from khmer._oxli.labeling import GraphLabels

from khmer._oxli.legacy_partitioning import SubsetPartition, PrePartitionInfo

from khmer._oxli.parsing import (FastxParser, SanitizedFastxParser,
                                 BrokenPairedReader)

from khmer._oxli.readaligner import ReadAligner

from khmer._oxli.utils import get_n_primes_near_x, is_prime, FILETYPES
from khmer._oxli.utils import get_version_cpp as __version_cpp__

import sys


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


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

