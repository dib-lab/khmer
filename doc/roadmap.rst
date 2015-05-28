.. vim: set filetype=rst

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


