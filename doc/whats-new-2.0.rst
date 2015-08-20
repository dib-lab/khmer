.. vim: set filetype=rst

************************
What's New In khmer 2.0?
************************

New behavior
============

Digital normalization script now supports mixed paired and unpaired read input
------------------------------------------------------------------------------

:program:`normalize-by-median.py` now supports mixed paired and unpaired (or
"broken-paired") input. Behavior can be forced to either treat all
reads as singletons or to require all reads be properly paired using
:option:`--force-single` or :option:`--paired`, respectively. If
:option:`--paired` is set, :option:`--unpaired-reads` can be used to
include a file of unpaired reads. The unpaired reads will be examined
after all of the other sequence files.

Reservoir sampling script extracts paired reads by default
----------------------------------------------------------

:program:`sample-reads-randomly.py` now retains pairs in the output, by
default.  This can be overridden to match previous behavior
with :option:`--force_single`.

New scripts
===========

Estimate number of unique kmers
-------------------------------

`unique-kmers.py` estimates the k-mer cardinality of a dataset using the
HyperLogLog probabilistic data structure. This allows very low memory
consumption, which can be configured through an expected error rate.
Even with low error rate (and higher memory consumption), it is still much
more efficient than exact counting and alternative methods.
It supports multicore processing (using OpenMP) and streaming,
and so can be used in conjunction with other scripts (like
`normalize-by-median.py` and `filter-abund.py`).

Incompatible changes
====================

New datastructure and script names
----------------------------------

For clarity the Count-Min Sketch based data structure previously known as
"counting_hash" or "counting_table" and variations of these is now known as
``countgraph``. Likewise with the Bloom Filter based data structure previously
known at "hashtable", "presence_table" and variations of these is now known as
``nodegraph``. Many options relating to 'table' have been changes to 'graph'.
Some scripts have been renamed: ``load-into-counting.py`` is now
:program:`load-into-counting.py`; ``load-graph.py`` is now
:program:`load-into-graph.py`.

New parameter for tablesize/number of table parameters.
-------------------------------------------------------

There is now a :option:`-M`/:option:`--max-memory-usage` parameter
that sets the number of tables (:option:`-N`/:option:`--num_tables`)
and tablesize (:option:`-x`/:option:`--max-tablesize`) parameters
automatically to match the desired memory usage.

(:option:`--min-tablesize` was also renamed to
:option:`--max-tablesize` to reflect this more desirable behavior.)

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

Scripts now output columnar data in CSV format by default
---------------------------------------------------------

All scripts that output any kind of columnar data now do so in CSV format,
with headers.  Previously this had to be enabled with :option:`--csv`.
(Affects :program:`abundance-dist-single.py`, :program:`abundance-dist.py`,
:program:`count-median.py`, and :program:`count-overlap.py`.)
:program:`normalize-by-median.py` also now outputs CSV when :option:`-R` is
used.

load-graph.py no longer appends .pt to the specified filename
-------------------------------------------------------------

Previously, `load-graph.py` appended a `.pt` extension to the
specified output filename and partition-graph appended a `.pt` to the
given input filename.  Now, `load-graph.py` writes to the specified
output filename and `partition-graph.py` does not append a `.pt` to
the given input filename.

Removed script
--------------

``count-overlap.py`` has been removed.
