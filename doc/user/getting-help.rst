..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2015 Michigan State University
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

===============
How to get help
===============

First, be sure that you:

#. Read the documentation (this site)

#. Google search for the error output and/or keywords related to your problem.
   Here you can search results from the mailing list, where others may
   have discussed solutions to the same issue.

.. raw:: html

    <form action="http://google.com/search" method="get">
    <input type="text" name="q" size="28" maxlength="255" value="" />
    <input type="submit" value="Google Search" />
    <br/>
    <input type="checkbox" name="sitesearch"
    value="http://lists.idyll.org/pipermail/khmer" checked /> only search
    khmer discussion email archive<br/>
    </form>

Mailing List
------------

The primary way to get help is through the khmer discussion list:
http://lists.idyll.org/listinfo/khmer, though we are also available for
closer-to-realtime support via `Gitter <https://gitter.im/dib-lab/khmer>`_. 

Asking a question
-----------------

#. Include your:

   * OS version (Mac OS X or Linux):  ``uname -mrs``
   * Python version:  ``python --version``
   * and khmer version:  ``pip freeze | grep khmer``

#. Precisely describe what you are trying to do.  Reread it from the
   perspective of someone else trying to reproduce your task.

#. Copy-and-paste the exact command that is causing the problem.  Include the
   steps you performed leading up to the issue.

#. Include the complete error message; if it is large, include a link to a
   file.

GitHub
------

You are also welcome to report an issue you are having using GitHub:
https://github.com/dib-lab/khmer/issues/new
