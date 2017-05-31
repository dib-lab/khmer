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

   Contact: khmer-project@idyll.org

A quick guide to testing (for khmer)
====================================

This document is for contributors new to automated testing, and explains
some of the motivation and logic behind the khmer project's testing
approach.

----

One of our most important "secret sauces" for khmer development is
that we do a fair bit of testing to make sure our code works and keeps
working!

* We maintain fairly complete test coverage of our code.  What this
  means is that we have automated tests that, when run, execute most
  of the lines of Python and C++ code in our src/, khmer/ and scripts/
  directories.  This doesn't *guarantee* things are correct, but it
  does mean that at least most of the code works at some basic level.

* we have other tests that we run periodically (for example, before
  each release) -- see :doc:`release` for details.  These tests
  check that our code works on multiple systems and with other
  people's software.

CTB and others have written a great deal about testing, and testing in
Python in particular.  Here's an `introductory guide
<http://ivory.idyll.org/articles/nose-intro.html>`__ CTB wrote a long
time ago.  You might also be interested in reading `this description
of the different kinds of tests
<http://www.ibm.com/developerworks/library/j-test/index.html>`__.

For the more general motivation, see `the Lack of Testing Death Spiral
<http://ivory.idyll.org/blog/software-quality-death-spiral.html>`__.

But... how do you do testing??

----

First, let's talk about specific goals for testing.  What should you
be aiming for tests to do?  You can always add more testing code, but
that might not be useful if they are redundant or over-complicated.

An overall rule is to "keep it simple" -- keep things as simple as
possible, testing as few things as possible in each test.

We suggest the following approach to writing tests for **new code**:

#. Write a test that just *runs* the new code, generally by copying existing
   test code to a new test and changing it.  Don't do anything clever for the
   first test -- just run something straightforward, and try to use existing
   data.

#. Decide which use cases should be tested.  This is necessarily code
   specific but our main advice is "don't be clever" -- write some tests
   to make sure that the code basically works.

#. Add in tests for edge cases.  By this we mean look for special cases in
   your code -- if statements, fence-post bound errors, etc. -- and write
   tests that exercise those bits of code specifically.

#. Make sure tests that expect a function call to fail (esp. with
   fail_ok=True) are failing for the expected reason. Run the code from the
   command line and see what the behavior is. For troubleshooting tests,
   catch the error with try: ... except: or print err.

For adding tests to **old code**, we recommend a mix of two approaches:

#. use `"stupidity driven testing"
   <http://ivory.idyll.org/blog/stupidity-driven-testing.html>`__ and
   write tests that recapitulate bugs before we fix those bugs.

#. look at test coverage (see `khmer's cobertura test coverage, here
   <http://ci.ged.msu.edu/job/khmer-master/label=linux/cobertura>`__) and
   identify lines of C++ or Python code that are not being executed by
   the current tests.  Then write new tests targeting the new code.

----

Next, to add a test, you have two options: either write a new one from
scratch, or copy an existing one.  (We recommend the latter.)

To write a new one, you'll need to know how to write tests. For
getting an idea of the syntax, read this `introductory guide
<http://ivory.idyll.org/articles/nose-intro.html>`__ and the `writing tests
documentation from Astropy
<http://docs.astropy.org/en/v1.1/development/testguide.html#writing-tests>`__.
Then find the right file in ``tests/*.py`` and add your test!

A better approach is, frankly, to go into the existing test code, find
a test that does something similar to what you want to do, copy it,
rename it, and then modify it to do the new test.

----

Finally, *where* do you add new tests and how do you run just *your* test?

Put new tests somewhere in ``tests/*.py``.  If you have trouble
figuring out what file to add them to, just put them in *some* file
and we'll help you figure out where to move them when we do code
review.

To run one specific test rather than all of them, you can do::

  py.test tests/test_scripts.py::test_load_into_counting

Here, you're running just one test -- the test function named
``test_load_into_counting`` in the file ``test_scripts.py``.

You can also invoke the test via setup.py, which is a bit more verbose::

  ./setup.py test --addopts "tests/test_scripts.py::test_load_into_counting"

----

Let's consider a simple test as an example.
The following code ensures that a k-mer and its reverse complement hash to the same value, since they represent the same molecule (just observed from a different orientation).

.. code:: python

    def test_kmer_rc_same_hash():
        kmer = 'GATTACAGATTACAGATTACA'
        kmer_rc = 'TGTAATCTGTAATCTGTAATC'

        ct = Counttable(21, 1e5, 2)
        assert ct.hash(kmer) == ct.hash(kmer_rc)

----

This example tests only a single function.
Tests that execute entire scripts and tests involving file I/O can require a bit more code.
Fortunately, we've created some helper functions that make this quite a bit easier in khmer.
See ``tests/test_scripts.py`` for some examples of code for executing scripts, capturing their output, and tidying up afterwards.
Also, see ``tests/khmer_tst_utils.py`` for helper functions that assist with loading test data and creating temporary output files.
