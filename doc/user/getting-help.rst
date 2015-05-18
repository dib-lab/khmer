.. vim: set filetype=rst

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
http://lists.idyll.org/listinfo/khmer

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

You are also welcome to report an issue you are having using GitHub::
https://github.com/dib-lab/khmer/issues/new
