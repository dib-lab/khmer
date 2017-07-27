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


khmer/Oxli Binary File Formats
==============================

- C++ macro definitions are given in parenthesis.
- C++ types are given in square brackets.
- ``Len`` is the field's size, in bytes, and ``Off`` is the field's zero-based
  byte offset in the file/section.

khmer v1.4 and previous
~~~~~~~~~~~~~~~~~~~~~~~

CountingHash
------------

The header is in the format below, in file offset order. There is no magic
string.

================== =========== ==============================================
Field               Length      Value
================== =========== ==============================================
Version             1           ``0x04`` (``SAVED_FORMAT_VERSION``)
File Type           1           ``0x01`` (``SAVED_COUNTING_HT``)
Use Bigcount        1           ``1`` if bigcounts is used, else ``0``
K-size              1           k-mer length, ``1 <= k <= 32``
Number of Tables    1           Number of Count-min Sketch tables
================== =========== ==============================================


khmer v2.0 formats
~~~~~~~~~~~~~~~~~~


Magic string
------------

All formats shall have the "magic string" ``OXLI`` as their first bytes, after
any external compression/encoding (e.g. gzip encapsulation) is removed. Note
that this makes them incompatible with older versions of khmer.

Countgraph
----------

(a.k.a ``CountingHash``, a Count-min Sketch)

The header is in the format below, again in the order of file offset.

================== ===== ===== ==============================================
Field               Len   Off     Value
================== ===== ===== ==============================================
Magic string        4       0   ``OXLI`` (``SAVED_SIGNATURE``)
Version             1       4   ``0x04`` (``SAVED_FORMAT_VERSION``)
File Type           1       5   ``0x01`` (``SAVED_COUNTING_HT``)
Use Bigcount        1       6   ``0x01`` if bigcounts is used, else ``0x00``
K-size              4       7   k-mer length, ``ht._ksize``. [``uint32_t``]
Number of Tables    1      11   Number of Count-min Sketch tables,
                                ``ht._n_tables``. [``uint8_t``]
Occupied Bins       8      12   Number of occupied bins
================== ===== ===== ==============================================

Then follows the Countgraph's tables. For each table:

================== ===== ===== ==============================================
Field               Len   Off     Value
================== ===== ===== ==============================================
Table size          8       0   Length of this table, ``ht._tablesizes[i]``.
                                [``uint64_t``]
Bins                N       8   This table's bins, length given by previous
                                field. [``uint8_t``]
================== ===== ===== ==============================================

Then follows a single value, the [``uint64_t``] number of ``kmer: count``
pairs. Then follows the Bigcount map, if this number is greater than zero. For
each kmer:

================== ===== ===== ==============================================
Field               Len   Off     Value
================== ===== ===== ==============================================
Kmer                8       0   Kmer's hash [``HashIntoType/uint64_t``].
Count               2       8   Kmer's count [``uint16_t``].
================== ===== ===== ==============================================


Nodegraph
---------

(a.k.a ``HashBits``, a Bloom Filter)

The header is in the format below, again in the order of file offset. Value
macro definitions are given in parenthesis

================== ===== ===== ==============================================
Field               Len   Off     Value
================== ===== ===== ==============================================
Magic string        4       0   ``OXLI`` (``SAVED_SIGNATURE``)
Version             1       4   ``0x04`` (``SAVED_FORMAT_VERSION``)
File Type           1       5   ``0x02`` (``SAVED_HASHBITS``)
K-size              4       6   k-mer length, ``ht._ksize``. [``unsigned int``]
Number of Tables    1      10   Number of Nodegraph tables. ``ht._n_tables``.
                                [``uint8_t``]
Occupied Bins       8      11   Number of occupied bins
================== ===== ===== ==============================================

Then follows the Nodegraph's tables. For each table:

================== ======= ===== ==============================================
Field               Len     Off     Value
================== ======= ===== ==============================================
Table size          8         0   Length of table, **in bits** (``uint64_t``).
Bins                N/8+1     8   This table's bytes, length given by previous
                                  field, divided by 8, plus 1 (``uint8_t``).
================== ======= ===== ==============================================
