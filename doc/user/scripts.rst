..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2010-2015 Michigan State University
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

******************************
khmer's command-line interface
******************************

The simplest way to use khmer's functionality is through the command
line scripts, located in the `scripts/
<https://github.com/dib-lab/khmer/tree/stable/scripts>`__ directory of the
khmer distribution.  Below is our documentation for these scripts.  Note
that all scripts can be given ``-h``/``--help`` which will print out
a list of arguments taken by that script.

Scripts that use k-mer counting tables or k-mer graphs take an
:option:`-M <load-into-counting.py -M>` parameter, which sets the maximum
memory usage in bytes. This should generally be set as high as possible; see
:doc:`choosing-table-sizes` for more information.

1. :ref:`scripts-counting`
2. :ref:`scripts-partitioning`
3. :ref:`scripts-diginorm`
4. :ref:`scripts-read-handling`

.. note::
 
   Almost all scripts take in either FASTA and FASTQ format, and
   output the same.

   Gzip and bzip2 compressed files are detected automatically. 

.. _scripts-counting:

k-mer counting and abundance filtering
======================================

.. autoprogram:: load-into-counting:get_parser()
        :prog: load-into-counting.py

.. autoprogram:: abundance-dist:get_parser()
        :prog: abundance-dist.py

.. autoprogram:: abundance-dist-single:get_parser()
        :prog: abundance-dist-single.py

.. autoprogram:: filter-abund:get_parser()
        :prog: filter-abund.py

.. autoprogram:: filter-abund-single:get_parser()
        :prog: filter-abund-single.py

.. autoprogram:: trim-low-abund:get_parser()
        :prog: trim-low-abund.py

.. autoprogram:: count-median:get_parser()
        :prog: count-median.py

.. autoprogram:: unique-kmers:get_parser()
        :prog: unique-kmers.py

.. _scripts-partitioning:

Partitioning
============

.. autoprogram:: do-partition:get_parser()
        :prog: do-partition.py

.. autoprogram:: load-graph:get_parser()
        :prog: load-graph.py

See :program:`extract-partitions.py` for a complete workflow.

.. autoprogram:: partition-graph:get_parser()
        :prog: partition-graph.py

See 'Artifact removal' to understand the stoptags argument.

.. autoprogram:: merge-partitions:get_parser()
        :prog: merge-partition.py

.. autoprogram:: annotate-partitions:get_parser()
        :prog: annotate-partitions.py

.. autoprogram:: extract-partitions:get_parser()
        :prog: extract-partitions.py
 
Artifact removal
----------------

The following scripts are specialized scripts for finding and removing
highly-connected k-mers (HCKs).  See :doc:`partitioning-big-data`.

.. autoprogram:: make-initial-stoptags:get_parser()
        :prog: make-initial-stoptags.py

.. autoprogram:: find-knots:get_parser()
        :prog: find-knots.py

.. autoprogram:: filter-stoptags:get_parser()
        :prog: filter-stoptags.py

.. _scripts-diginorm:

Digital normalization
=====================

.. autoprogram:: normalize-by-median:get_parser()
        :prog: normalize-by-median.py

.. _scripts-read-handling:

Read handling: interleaving, splitting, etc.
============================================

.. autoprogram:: extract-long-sequences:get_parser()
        :prog: extract-long-sequences.py

.. autoprogram:: extract-paired-reads:get_parser()
        :prog: extract-paired-reads.py

.. autoprogram:: fastq-to-fasta:get_parser()
        :prog: fastq-to-fasta.py

.. autoprogram:: interleave-reads:get_parser()
        :prog: interleave-reads.py

.. autoprogram:: readstats:get_parser()
        :prog: readstats.py

.. autoprogram:: sample-reads-randomly:get_parser()
        :prog: sample-reads-randomly.py

.. autoprogram:: split-paired-reads:get_parser()
        :prog: split-paired-reads.py
