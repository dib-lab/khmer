A quick guide to the khmer codebase
===================================

This document describes the khmer project's directory layout.

----

The ChangeLog file lists changes to the codebase, most recent first.

The lib/ directory contains all of the C++ code.

The khmer/ directory contains the khmer package (khmer/__init__.py etc)
and the C++-to-Python bridge (khmer/_khmermodule.cc).

The scripts/ and sandbox/ directory contain Python command-line scripts.

The tests/ directory contains all of the tests.  Each test is a function in
one of the tests/test*.py files.
