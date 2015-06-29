.. vim: set filetype=rst

What's New In khmer 2.0?
########################

Incompatible changes
====================

New parameter for tablesize/number of table parameters.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is now a :option:`-M`/:option:`--max-memory-usage` parameter
that sets the number of tables (:option:`-N`/:option:`--num_tables`)
and tablesize (:option:`-x`/:option:`--max-tablesize`) parameters
automatically to match the desired memory usage.

(:option:`--min-tablesize` was also renamed to
:option:`--max-tablesize` to reflect this more desirable behavior.)

Binary file formats have changed!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All binary khmer formats (presence tables, counting tables, tag sets,
stop tags, and partition subsets) have changed. Files are now
pre-pended with the string ``OXLI`` to indicate that they are from
this project.

Files of the above types made in previous versions of khmer are not compatible
with v2.0; the reverse is also true.
