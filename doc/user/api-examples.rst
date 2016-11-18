..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2016 The Regents of the University of California.
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

Example API Usage
=================

Examples of how to run khmer command-line scripts are included in the main khmer documentation, as well as the `khmer protocols collection <http://khmer-protocols.readthedocs.io>`__.
This document is intended for programmers and developers who want to use khmer's Python or C++ API.

**Note**: The khmer command-line scripts are under `semantic versioning <http://semver.org/>`__, but (as of khmer v2) neither the Python API nor the C++ API are covered under semantic versioning.
Consequently, updates to the khmer code base may require changes to the examples below at any time.
We expect to have the Python API under semantic versioning by version 5.0.


Python API
----------

The following example demonstrates the basics of counting k-mers with khmer.

.. code:: python

    >>> import khmer
    >>>
    >>> ksize = 11
    >>>
    >>> # exact counting with and odd k requires a single table with 4^k bytes
    >>> # even k will collapse k-mers with their reverse complement, and requires 4^(k-1) + k bytes
    >>> cg = khmer.Countgraph(ksize, 4**ksize, 1)
    >>>
    >>> # load all k-mers from the given string (returns the number of k-mers processed)
    >>> cg.consume('ATGGCGATGGCAAGTAGGACCCAGATGGACCAAAG')
    25
    >>>
    >>> # what is the count of the k-mer "ATGGCGATGGC"?
    >>> cg.get('ATGGCGATGGC')
    1
    >>> # increase the count of the k-mer by 1
    >>> cg.add('ATGGCGATGGC')
    1
    >>> # what is the k-mer's count now?
    >>> cg.get('ATGGCGATGGC')
    2
    >>>
    >>> # what is the count of this bogus k-mer?
    >>> cg.get('GTGGCGATGGC')
    0
    >>>

C++ API
-------

The `count-stuff.cpp` file in the `examples/c++-api/` directory contains source code that calls khmer's C++ API directly.
To run this example, issue the following commands from your terminal.
The `Makefile` in `examples/c++-api/` shows how to compile and link against the khmer library.

.. code:: bash

    cd examples/c++-api/
    make
    ./count-stuff
