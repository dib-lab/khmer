..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2014-2015 Michigan State University
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

***************************
Roadmap to v2.0, v3.0, v4.0
***************************

Background
==========

To make the khmer project easier to use and easier to build upon several
fundamental changes need to happen. This document outlines our plan to do so
while minimizing the impact of these changes on our existing users.

The discussion that lead to this document can be read at
https://github.com/dib-lab/khmer/issues/389

Remainder of v1.x series
========================

Start of transition to a single entrypoint named `oxli`. This will be exempt
from the project's semantic versioning and will be advertised as experimental
and unstable.

Migration of Python script functionality to a Python module named `oxli`. As
the code moves over there will be no change to external script functionality or
their command line interfaces.

v2.x series
===========

`oxli` command is now under semantic versioning. Scripts are still the
advertised and preferred entry point for users. Developers and workflow systems
can start to trial `oxli` but need not switch until 3.0. New functionality is
added to both the scripts and the `oxli` command.

v3.0 and project renaming
=========================

Project renamed to 'oxli'; all references to 'khmer' removed from the code and
documentation except for a single note in the docs. All scripts dropped as
their functionality has been moved to the `oxli` command. Websites that we
maintain that have 'khmer' in the URL will have redirects installed.

Refinement of the Python API continues, however it is not part of the semantic
versioning of the project.

v4.0
====

The semantic versioning now extends to the Python API.

Python API wishlist
===================

Python 3.0 support

API for multiple container types and implementation of the same.

Cleanup of Python/C++ class hierarchy to cut down on boilerplate glue code.

Switch to new-style Python objects (see LabelHash & Hashbits)


