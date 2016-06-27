..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2014-2015 Michigan State University
   Copyright (C) 2015-2016 The Regents of the University of California.
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

Coding guidelines and code review checklist
===========================================

This document is for anyone who want to contribute code to the khmer
project, and describes our coding standards and code review checklist.

C++ standards
-------------

Any feature in C++11 is fine to use. Specifically we support features found in
GCC 4.8.2. See https://github.com/dib-lab/khmer/issues/598 for an in-depth
discussion.

Coding standards
----------------

All plain-text files should have line widths of 80 characters or less unless
that is not supported for the particular file format.

For C++, we use `Todd Hoff's coding standard
<http://www.possibility.com/Cpp/CppCodingStandard.html>`__, and
`astyle -A10 / "One True Brace Style"
<http://astyle.sourceforge.net/astyle.html>`__ indentation and
bracing.  Note: @CTB needs Emacs settings that work for this.

Vim users may want to set the ARTISTIC_STYLE_OPTIONS shell variable to "-A10
--max-code-length=80" and run ```:%!astyle``` to reformat. The four space
indentation can be set with::

	set expandtab
	set shiftwidth=4
	set softtabstop=4

For Python, `PEP 8 <http://www.python.org/dev/peps/pep-0008/>`__ is our
standard. The ```pep8``` and ```autopep8``` Makefile targets are helpful.

Code, scripts, and documentation must have their spelling checked.

Python-based `codespell` can be applied to multiple files easily. `codespell`
can be installed via the following::

        mkdir ~/bin
        git clone git@github.com:lucasdemarchi/codespell.git
        cd codespell
        make prefix=${HOME} install
        export PATH=$PATH:~/bin/

Note, if you want codespell to always be available you will need to add the
`export` line to your `${HOME}\.bashrc` or equivalent.

To run codespell over only what has been changed on the branch `my-branch`::

        git diff master..my-branch > diff_file
        codespell diff_file

To run codespell over a single file::

        codespell path/to/file

To make codespell fix the issues it finds automatically::

        codespell -w path/to/file

Please note that as `codespell` works off of a listing of possible
misspellings it may not catch all errors. If you find a spelling error that
is not caught by `codespell` feel free to open a pull request at the `project
page <https://github.com/lucasdemarchi/codespell>`_ to add it to the
dictionary.

Vim users can run::

        :setlocal spell spelllang=en_us

Use `]s` and `[s` to navigate between misspellings and `z=` to suggest a
correctly spelled word. `zg` will add a word as a good word.

GNU `aspell` can also be used to check the spelling in a single file::

        aspell check --mode ccpp $filename

Code Review
-----------

Please read `11 Best Practices for Peer Code Review
<http://smartbear.com/SmartBear/media/pdfs/WP-CC-11-Best-Practices-of-Peer-Code-Review.pdf>`__.

See also `Code reviews: the lab meeting for code
<http://fperez.org/py4science/code_reviews.html>`__ and
`the PyCogent coding guidelines
<http://pycogent.org/coding_guidelines.html>`__.

Checklist
---------

Each pull request should be automatically populated with the following
checklist::

   - [ ] Is it mergeable?
   - [ ] `make test` Did it pass the tests?
   - [ ] `make clean diff-cover` If it introduces new functionality in
     `scripts/` is it tested?
   - [ ] `make format diff_pylint_report cppcheck doc pydocstyle` Is it well
     formatted?
   - [ ] Did it change the command-line interface? Only additions are allowed
     without a major version increment. Changing file formats also requires a
     major version number increment.
   - [ ] Is it documented in the `ChangeLog`?
     http://en.wikipedia.org/wiki/Changelog#Format
   - [ ] Was a spellchecker run on the source code and documentation after
     changes were made?
   - [ ] Do the changes respect streaming IO? (Are they
     tested for streaming IO?)
   - [ ] Is the Copyright year up to date?

**Note** that after you submit the pull request you can check and uncheck
the individual boxes on the formatted comment; no need to put x or y
in the middle.

CPython Checklist
-----------------

Here's a checklist for new CPython types with future-proofing for Python 3::

   - [ ] the CPython object name is of the form `khmer_${OBJECTNAME}_Object`
   - [ ] Named struct with `PyObject_HEAD` macro
   - [ ] `static PyTypeObject khmer_${OBJECTNAME}_Type` with the following
     entries
      - [ ] `PyVarObject_HEAD_INIT(NULL, 0)` as the object init (this includes
        the `ob_size` field).
      - [ ] all fields should have their name in a comment for readability
      - [ ] The `tp_name` filed is a dotted name with both the module name and
        the name of the type within the module. Example: `khmer.ReadAligner`
      - [ ] Deallocator defined and cast to `(destructor)` in tp_dealloc
        - [ ] The object's deallocator must be
          `Py_TYPE(obj)->tp_free((PyObject*)obj);`
      - [ ] Do _not_ define a `tp_getattr`
      - [ ] BONUS: write methods to present the state of the object via
        `tp_str` & `tp_repr`
      - [ ] _Do_ pass in the array of methods in `tp_methods`
      - [ ] _Do_ define a new method in `tp_new`
   - [ ] PyMethodDef arrays contain doc strings
      - [ ] Methods are cast to `PyCFunctions`s
   - [ ] Type methods use their type Object in the method signature.
   - [ ] Type creation method decrements the reference to self
     (`Py_DECREF(self);`) before each error-path exit (`return NULL;`)
   - [ ] No factory methods. Example: `khmer_new_readaligner`
   - [ ] Type object is passed to `PyType_Ready` and its return code is checked
     in `MOD_INIT()`
   - [ ] The reference count for the type object is incremented before adding
     it to the module: `Py_INCREF(&khmer_${OBJECTNAME}_Type);`.
