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
However, the intrepid user may be interested in having more direct access to khmer via its Python or C++ API.

The ``examples/python-api`` and ``examples/c++-api`` directories contain several small programs that demonstrate how one would write a program that uses the khmer API.
These programs are run frequently as part of our continuous integration build, so they are quite stable.
Note, however, that unlike the khmer command-line scripts, the Python and C++ API **ARE NOT** under `semantic versioning <http://semver.org/>`__.
Consequently, internal changes to the khmer codebase may require changes to the API examples at any time.
We expect to have the Python API under semantic versioning :doc:`by version 5.0 <../roadmap>`.
