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

Guidelines for continued development
====================================

The khmer :doc:`Getting Started <getting-started>` documentation is intended for
first-time contributors, and is designed to make that first contribution as easy
and as painless as possible. For those with an interest in making continued
contributions (or those with an obligation to maintain and contribute to the
codebase), this document describes the coding guidelines we follow, as well as
some tips that will improve the development process for everyone involved.

Beyond your first contribution
------------------------------

If your :doc:`first contribution to khmer <getting-started>` has been
accepted and merged into the codebase, you may be wondering what to do next. You
can poke around at additional "low-hanging fruit" issues, but if you'd like to
sink your teeth into something more substantial, here are a few suggestions.

* If you're knowledgeable in C++ and/or Python and/or documentation
  and/or biology, we'd love to attract further contributions to khmer.
  Please visit the issues list and browse about and find something
  interesting looking.

* One general thing we'd like to do is increase our test coverage.
  You can go find test coverage information `at Codecov.io
  <https://codecov.io/gh/dib-lab/khmer>`__ by scrolling to the bottom and
  clicking on individual files; or, ask us on khmer-project@idyll.org for
  suggestions.

* Ask us! Ask khmer-project@idyll.org for suggestions on what to do next.
  We can suggest particularly ripe low-hanging fruit, or find some other
  issues that suit your interests and background.

* You can also help other people out by watching for new issues or
  looking at pull requests. Remember to be nice and polite!

Programming languages
---------------------

All Python code in khmer must run correctly in both Python version 2 and 3.

For C++ code, any feature in C++11 is fine to use. Specifically we support
features found in GCC 4.8.2. Our automated tests use gcc 4.8.4 on linux. See
https://github.com/dib-lab/khmer/issues/598 for an in-depth discussion. Please
do not use features from C++14 or newer.

Code style standards
--------------------

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

        pip install codespell

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

Cython Style
~~~~~~~~~~~~

Cython code can become messy very quickly, and as such, we have guidelines
for style and structure.

When wrapping code from liboxli:

- `extern` definition should begin with `Cp`; for example, `CpHashtable` wraps
  `oxli::Hashtable`.
- If the extension class wrapping the liboxli class stores it as a pointer,
  it should be named `_this`. If it wraps a stack object directly, it should
  be named `_obj`.

For imports,

- `libc` cimports next,
- then `libcpp` imports and cimports.
- followed by cimports
- and finally, regular imports.

Generally,

- Pure C methods should be underscore prefixed.
- `cpdef` methods do not need to be underscore prefixed.


Resolving merge conflicts
-------------------------

It is possible that when you do a `git pull` you will get a "merge
conflict" -- This is what happens when something changed in the branch you're
pulling in in the same place you made a change in your local copy. This
frequently happens in the `ChangeLog` file.

Git will complain loudly about merges and tell you specifically in which
files they occurred. If you open the file, you'll see something vaguely
like this in the place where the merge occurred::

   <<<<<<< HEAD
   Changes made on the branch that is being merged into. In most cases,
   this is the branch that you have currently checked out
   =======
   Changes made on the branch that is being merged in, almost certainly
   master.
   >>>>>>> abcde1234

Though there are a variety of tools to assist with resolving merge
conflicts they can be quite complicated at first glance and it is usually
easy enough to manually resolve the conflict.

To resolve the conflict you simply have to manually 'meld' the changes
together and remove the merge markers.

After this you'll have to add and commit the merge just like any other set
of changes. It's also recommended that you run tests.


Virtual environments
--------------------

The khmer package, like many software packages, relies on other third-party
software. Some of this software has been bundled together with khmer and is
compiled when you invoke ``make`` on the command line. But some of the software
khmer depends on is distributed as Python packages separately from khmer.

Python `virtual environments <https://pypi.python.org/pypi/virtualenv>`_ were
designed to isolate a stable development environment for a particular project.
This makes it possible to maintain different versions of a Python package for
different projects on your computer.

The installation instructions in the :doc:`Getting Started <getting-started>`
docs install the ``virtualenv`` command on your computer. After completing those
instructions, you can create a virtual environment with the command::

    virtualenv -p python2 env/

(You can substitute `python3` for `python2` if Python version 3 is installed on
your system.) This command will create a new directory `env/` containing your
new virtual environment. The command::

    source env/bin/activate

will activate the virtual environment. Now any Python packages that you install
with ``pip`` or ``make install-dep`` will be installed into your isolated
virtual environment.

Note that any time you create a new terminal session, using the virtual
environment requires that you re-activate it.

Pull request cleanup (commit squashing)
---------------------------------------

Submitters are invited to reduce the numbers of commits in their pull requests
either via `git rebase -i dib/master` or this recipe::

        git pull # make sure the local is up to date
        git pull dib master # get up to date
        # fix any merge conflicts
        git status # sanity check
        git diff dib/master # does the diff look correct? (no merge markers)
        git reset --soft dib/master # un-commit the differences from dib/master
        git status # sanity check
        git commit --all # package all differences in one commit
        git status # sanity check
        git push # should fail
        git push --force # override what's in GitHub's copy of the branch/pull request


Code Review
-----------

Please read `11 Best Practices for Peer Code Review
<http://smartbear.com/SmartBear/media/pdfs/WP-CC-11-Best-Practices-of-Peer-Code-Review.pdf>`__.

See also `Code reviews: the lab meeting for code
<http://fperez.org/py4science/code_reviews.html>`__ and
`the PyCogent coding guidelines
<http://pycogent.org/coding_guidelines.html>`__.

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
