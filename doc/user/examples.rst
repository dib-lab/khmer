..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2013-2015 Michigan State University
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

**************
A few examples
**************

See the `examples <https://github.com/dib-lab/khmer/tree/stable/examples>`__
directory for complete examples.

STAMPS data set
===============

The 'stamps' data set is a fake metagenome-like data set containing
two species, mixed at a 10:1 ratio.  The source genomes are
in `data/stamps-genomes.fa
<https://github.com/dib-lab/khmer/tree/stable/data/stamps-genomes.fa>`__. 
The reads file is in `data/stamps-reads.fa.gz
<https://github.com/dib-lab/khmer/tree/stable/data/stamps-reads.fa.gz>`__,
and consists of 100-base reads with a 1% error rate.

The example shows how to construct k-mer abundance histograms, as well
as the effect of digital normalization and partitioning on the k-mer
abundance distribution.

See `the script for running everything
<https://github.com/dib-lab/khmer/blob/stable/examples/stamps/do.sh>`__
and `the IPython Notebook
<http://nbviewer.ipython.org/urls/raw.github.com/dib-lab/khmer/stable/examples/stamps%2520k-mer%2520distributions.ipynb>`__.

For an overall discussion and some slides to explain what's going on,
visit `the Web site for a 2013 HMP metagenome assembly webinar that
Titus Brown gave <http://ged.msu.edu/angus/2013-hmp-assembly-webinar/exploring-stamps-data.html>`__.

