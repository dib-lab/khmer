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

Deploying the khmer project tools on Galaxy
===========================================

This document is for people interested in deploying the khmer tools on
the Galaxy platform.

----

We are developing the support for running all the khmer scripts in `Galaxy
<http://galaxyproject.org/>`__.

Install the tools & tool description
------------------------------------

In the administrative interface select "Search and browse tool sheds" under
the heading "Tool sheds". Click on "Galaxy test tool shed" and search for
khmer. Click on the "khmer" button and choose "Preview and install". Click the
"Install to Galaxy" button at the top. At the bottom of the next page click
the "Install" button.

Single Output Usage
-------------------

For one or more files into a single file:

#. Choose 'Normalize By Median' from the 'khmer protocols' section of the
   'Tools' menu.

#. Compatible files already uploaded to your Galaxy instance should be listed.
   If not then you may need to `set their datatype manually
   <https://wiki.galaxyproject.org/Learn/Datatypes>`__.

#. After selecting the input files specify if they are paired-interleaved
   or not.

#. Specify the sample type or show the advanced parameters to set the tablesize
   yourself. Consult :doc:`../user/choosing-table-sizes` for assistance.
