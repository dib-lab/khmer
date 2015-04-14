.. vim: set filetype=rst

===============
How to get help
===============

**First**, be sure that you:

#. Read the documentation (this site)

#. Work through relevant examples in the `khmer-recipes repository <https://github.com/ged-lab/khmer-recipes>`__. Find **khmer-recipes** `documentation here <http://khmer-recipes.readthedocs.org/en/latest/#>`__.

#. Google search for the error output and/or keywords related to your problem.  This will often bring up results from the mailing list, where others may have discussed solutions to the same issue.

.. raw:: html

    <form action="http://google.com/search" method="get">
    <input type="text" name="q" size="28" maxlength="255" value="" /> <input type="submit" value="Google Search" />
    <br/>
    <input type="checkbox" name="sitesearch" value="http://lists.idyll.org/pipermail/khmer" checked /> only search khmer discussion email archive<br/>
    </form>

Mailing List
------------

The primary way to get help is through the khmer discussion list:
http://lists.idyll.org/listinfo/khmer

When you ask a question, it is important that you ask in a way that ensures you get the best answer possible and minimizes back-and-forth with the good people trying to help.

Asking a good question
----------------------

#. Be very specific in your subject line.
  - If well written, it will more likely catch the eye of someone who knows the solution.
  - Subject lines like *"arrgh cant get it to work"* or *"NEED HELP, POSTER DUE TOMORROW!!!!!"* are not likely to inspire anyone to help.
  - Example of a good subject line: *"Error message from normalize-by-median 'Hash writing file access failure:: Bad address'"*

2. Include your:

 a. OS version:  ``uname -mrs``

 b. Python version:  ``python --version``

 c. and khmer version:  ``pip freeze | grep khmer``

3. Precisely describe what you are trying to do.  Reread it from the perspective of someone else trying to reproduce your task.

#. Copy-and-paste the exact command that is causing the problem.  Include the steps you performed leading up to the issue.

#. Include the complete error message; if it is large, include a link to a file.

GitHub
------

If you think you've found a legitimate bug or other problem with
functionality, use GitHub to report an issue:
https://github.com/ged-lab/khmer/issues/new
