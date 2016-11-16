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

Command line scripts, ``scripts/``, and ``sandbox/``
====================================================

.. note::

   This document applies through khmer/oxli 2.0/3.0 (see
   :doc:`../roadmap`) - we will revisit when the Python API falls
   under semantic versioning for oxli 4.0.

khmer has two conflicting goals: first, we want to provide a reliable
piece of software to our users; and second, we want to be flexible and
enable exploration of new algorithms and programs.  To this end,
we've split our command line scripts across two directories,
``scripts/`` and ``sandbox/``.  The former is the staid, boring, reliable
code; the latter is a place for exploration.

As a result, we are committed to high test coverage, stringent code
review, and `Semantic Versioning <http://semver.org/>`__ for files in
``scripts/``, but explicitly *not* committed to this for files and
functionality implemented in ``sandbox/``.  So, putting a file into
``scripts/`` is a big deal, especially since it increases our maintenance
burden for the indefinite future.

We've roughed out the following process for moving scripts into ``scripts/``:

* Command line scripts start in ``sandbox/``;
* Once their utility is proven (in a paper, for example), we can propose to
  move them into ``scripts/``;
* There's a procedure for moving scripts from ``sandbox/`` into ``scripts/``.

Read on!

Sandbox script requirements and suggestions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All scripts in ``sandbox/`` must:

* be importable (enforced by ``test_import_all`` in
  ``test_sandbox_scripts.py``)
* be mentioned in ``sandbox/README.rst``
* have a hash-bang line (``#! /usr/bin/env python``) at the top
* be command-line executable (``chmod a+x``)
* have a Copyright message (see below)
* have lowercase names
* use '-' as a word separator, rather than '_' or CamelCase

All *new* scripts being added to ``sandbox/`` should:

* have decent automated tests
* be used in a protocol (see khmer-protocols) or a recipe (see khmer-recipes)
* be pep8 clean and pylint clean-ish (see ``make pep8`` and ``make_diff_pylint``).

Standard and reserved command line options for ``scripts/``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following options are reserved, and the short option flag cannot be
redefined in any official script.

* ``-h|--help`` - this must always print out a descriptive usage statement
* ``-v|--version`` - this must always print out the khmer version number
* ``-x|--max-tablesize`` - this must always specify the approximate table size
  for storing sketches of k-mer hashes
* ``-N|--n_tables`` - this must always specify the number of tables for storing
  sketches of k-mer hashes
* ``-M|--max-memory-usage`` - this must always specify the maximum amount of
  memory to be consumed for storing sketches of k-mer hashes
* ``-U|--unique-kmers`` - this must always specify the approximate number of
  unique k-mers in the data set
* ``-k|--ksize`` - this must always specify the k-mer size to use
* ``-q|--quiet`` - this must always indicate that the script's diagnostic output
  should be minimized or altogether eliminated

Additionally, all scripts in ``scripts/`` should have the following options.

* ``-h|--help``
* ``-v|--version``
* ``-f|--force`` - if applicable, override any sanity checks that may prevent
  the script from running

If an option is of type ``type=argparse.FileType('w')`` then you need to also
specify a ``metavar`` for the documentation and help formatting. Example::

    parser.add_argument('-R', '--report', metavar='report_filename',
        type=argparse.FileType('w'))

Copyright message
~~~~~~~~~~~~~~~~~

The copyright message should be of the form::

   #
   # This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   # Copyright (C) ___WHO___, ___YEAR(s)___. It is licensed under
   # the three-clause BSD license; see doc/LICENSE.txt.
   # Contact: khmer-project@idyll.org
   #

Where ___WHO___ is replaced with one or more of "Michigan State University" or
"The Regents of the University of California" and ___YEAR(s)___ is replaced
with the year or years the file was created or modified. The copyright
statement for new files should only refer to "The Regents of the University of
California".

Upgrading a script from 'sandbox' to 'scripts'
----------------------------------------------

First, everything needed (all library support code) should be already
committed to khmer master after the usual review process; the relevant
script(s) should be in ``sandbox/``.

Second, an issue should be started explicitly to discuss whether the
script(s) should be moved from ``sandbox/`` into ``scripts/``.  This issue
should discuss the general need for this script, outside of a particular
paper pipeline.  (Note that there is no imperative to move a script
out of ``sandbox/``; if we think it's useful code to have around and
want to keep it functioning, we should just add in automated tests and
otherwise level it up.)

Third, assuming we reach general agreement about moving the script(s)
into ``scripts/``, start a pull request to do so, referencing the
issue and containing the following checklist.  The PR should start by
moving the script from ``sandbox/`` into ``scripts/``, and moving the
tests out of the ``test_sandbox_scripts.py`` file.

Last but not least, intensive code review may raise more general
issues that could apply to the entire code base; if contentious or
needing discussion, these issues may be punted to general issues so as
to not block a merge.

A checklist for moving a script into the scripts/ directory from sandbox/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy or paste this checklist into the PR, in addition to the normal
development/PR checklist::

   - [ ] most or all lines of code are covered by automated tests (see output of ``make diff-cover``)
   - [ ] ``make diff_pylint`` is clean
   - [ ] the script has been updated with a ``get_parser()`` and added to doc/user/scripts.txt
   - [ ] argparse help text exists, with an epilog docstring, with examples and options
   - [ ] standard command line options are implemented
   - [ ] version and citation information is output to STDERR (`khmer_args.info(...)`)
   - [ ] support '-' (STDIN) as an input file, if appropriate
   - [ ] support designation of an output file (including STDOUT), if appropriate
   - [ ] script reads and writes sequences in compressed format
   - [ ] runtime diagnostic information (progress, etc.) is output to STDERR
   - [ ] script has been removed from sandbox/README.rst
