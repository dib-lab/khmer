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

Incompatible changes
====================

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
