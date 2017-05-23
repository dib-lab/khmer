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

*********************************
Contributors and Acknowledgements
*********************************

khmer is a product of the Lab for Data Intensive Biology at the University of
California, Davis (the succesor to the GED lab at Michigan State University),

   http://ivory.idyll.org/lab/

---

C. Titus Brown <titus@idyll.org> wrote the initial Bloom filter and
Count-Min Sketch implementations, has contributed to their continuing
maintenance and improvement and refactoring, and has contributed
extensively more generally to the codebase.

Jason Pell implemented many of the C++ k-mer filtering functions.

Qingpeng contributed code to do unique k-mer counting.

Adina Howe, Rosangela Canino-Koning, and Arend Hintze contributed
significantly to discussions of approaches and algorithms; Adina wrote
a number of scripts.

Jared T. Simpson (University of Cambridge, Sanger Institute) contributed
paired-end support for digital normalization.

Eric McDonald thoroughly revised many aspects of the code base, made
much of the codebase thread safe, and otherwise improved performance
dramatically.

Michael R. Crusoe took over maintainership June, 2013, and streamlined
and improved many of khmer's development and community processes.

Luiz Irber implemented an efficient HyperLogLog-based cardinality
estimator and has contributed extensively more generally to the
codebase.

Camille Scott has contributed significantly to khmer's assembly and
graph traversal functionality, and has contributed extensively more
generally to the codebase.

Tim Head has contributed to refactoring the core data structures,
performance benchmarking, and has contributed extensively more generally
to the codebase.

Daniel Standage took over maintainership May, 2016, has refined the
documentation extensively, contributed Python and C++ code examples,
refactored core data structures for a more extensible sequence loading
functionality, and contributed more generally to the codebase.

Last updated by DSS on 2017-05-22
