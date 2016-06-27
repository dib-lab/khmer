..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2015 The Regents of the University of California.
   It is licensed under the three-clause BSD license; see LICENSE.
   Contact: khmer-project@idyll.org

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Contact: khmer-project@idyll.org

************************
What's New In khmer 2.0?
************************

New behavior
============

Streaming I/O from Unix Pipes
-----------------------------

All scripts now accept input from named (like ``/dev/stdin``, or that created
using ``<( list )`` process substituion) and unamed pipes (like output piped in
from another program with ``|``). The STDIN stream can also be specified using
a single dash: ``-``.

New parameter for memory usage, and/or tablesize/number of table parameters.
----------------------------------------------------------------------------

There is now a :option:`-M <load-into-counting.py -M>`/
:option:`--max-memory-usage <load-into-counting.py --max-memory-usage>`
parameter that sets the number of tables (
:option:`-N <load-into-counting.py -N>`/
:option:`--n_tables <load-into-counting.py --n_tables>`) and tablesize
(:option:`-x <load-into-counting.py -x>`/:option:`--max-tablesize
<load-into-counting.py --max-tablesize>`) parameters automatically to match the
desired memory usage.

Digital normalization script now supports mixed paired and unpaired read input
------------------------------------------------------------------------------

:program:`normalize-by-median.py` now supports mixed paired and unpaired (or
"broken-paired") input. Behavior can be forced to either treat all
reads as singletons or to require all reads be properly paired using
:option:`--force_single <normalize-by-median.py --force_single>` or
:option:`--paired <normalize-by-median.py --paired>`, respectively. If
:option:`--paired <normalize-by-median.py --paired>` is set,
:option:`--unpaired-reads <normalize-by-median.py --unpaired-reads>` can be
used to include a file of unpaired reads. The unpaired reads will be examined
after all of the other sequence files.
:option:`normalize-by-median.py --quiet` can be used to reduce the amount of
diagnostic output.

Mixed-pair sequence file format support
---------------------------------------

:option:`split-paired-reads.py --output-orphaned`/:option:`-0
<split-paired-reads.py -0>` has been added to allow for orphaned reads and give
them a file to be sorted into.

Scripts now output columnar data in CSV format by default
---------------------------------------------------------

All scripts that output any kind of columnar data now do so in CSV format,
with headers.  Previously this had to be enabled with ``--csv``.
(Affects :program:`abundance-dist-single.py`, :program:`abundance-dist.py`,
:program:`count-median.py`, and :program:`count-overlap.py`.)
:option:`normalize-by-median.py --report` also now outputs in CSV format.

Reservoir sampling script extracts paired reads by default
----------------------------------------------------------

:program:`sample-reads-randomly.py` now retains pairs in the output, by
default.  This can be overridden to match previous behavior
with :option:`--force_single <sample-reads-randomly.py --force_single>`.

New scripts
===========

Estimate number of unique kmers
-------------------------------

:program:`unique-kmers.py` estimates the k-mer cardinality of a dataset using
the HyperLogLog probabilistic data structure. This allows very low memory
consumption, which can be configured through an expected error rate.
Even with low error rate (and higher memory consumption), it is still much
more efficient than exact counting and alternative methods.
It supports multicore processing (using OpenMP) and streaming,
and so can be used in conjunction with other scripts (like
:program:`normalize-by-median.py` and :program:`filter-abund.py`). This is the
work of Luiz Irber and it is the subject of a paper in draft.

Incompatible changes
====================

New datastructure and script names
----------------------------------

For clarity the Count-Min Sketch based data structure previously known as
"counting_hash" or "counting_table" and variations of these is now known as
``countgraph``. Likewise with the Bloom Filter based data structure previously
known at "hashbits", "presence_table" and variations of these is now known as
``nodegraph``. Many options relating to ``table`` have been changed to
``graph``.


Binary file formats have changed
--------------------------------

All binary khmer formats (presence tables, counting tables, tag sets,
stop tags, and partition subsets) have changed. Files are now
pre-pended with the string ``OXLI`` to indicate that they are from
this project.

Files of the above types made in previous versions of khmer are not compatible
with v2.0; the reverse is also true.

In addition to the ``OXLI`` string, the Nodegraph and Countgraph file format
now includes the number of occupied bins. See :doc:`dev/binary-file-formats`
for details.

load-graph.py no longer appends .pt to the specified filename
-------------------------------------------------------------

Previously, :program:`load-graph.py`` appended a ``.pt`` extension to the
specified output filename and :program:`partition-graph.py` appended a ``.pt``
to the given input filename.  Now, :program:`load-graph.py` writes to the
specified output filename and :program:`partition-graph.py` does not append a
``.pt`` to the given input filename.

Some reporting options have been turned always on
-------------------------------------------------

The total number of unique k-mers will always be reported every time a new
countgraph is made. The ``--report-total-kmers`` option has been removed from
:program:`abundance-dist-single.py`, :program:`filter-abund-single.py`, and
:program:`normalize-by-median.py` to reflect this. Likewise with
``write-fp-rate`` for :program:`load-into-counting.py` and
:program:`load-graph.py`; the false positive rate will always be
written to the ``.info`` files.

An uncommon error recovery routine was removed
----------------------------------------------

To simplify the codebase ``--save-on-failure`` and its helper option
``--dump-frequency`` have been removed from :program:`normalize-by-median.py`.

Single file output option names have been normalized
----------------------------------------------------

``--out`` is now ``--output`` for both :option:`normalize-by-median.py
<normalize-by-median.py --output>` and :option:`trim-low-abund.py
<trim-low-abund.py --output>`.

Miscellaneous changes
---------------------
The common option ``--min-tablesize`` was renamed to
:option:`--max-tablesize <load-into-counting.py --max-tablesize>` to reflect
this more desirable behavior.

In conjuction with the new :option:`split-paired-reads.py --output-orphaned`
option, the option ``--force-paired``/``-p`` has been eliminated.

As CSV format is now the default, the ``--csv`` option has been removed.

Removed script
--------------

`count-overlap.py
<http://khmer.readthedocs.org/en/v1.4.1/user/scripts.html#count-overlap-py>`__
has been removed.
