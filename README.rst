|Research software impact|
|Supported Python versions|
|khmer build status|
|Test coverage|
|BSD-3 licensed|

khmer
=====

Welcome to khmer: k-mer counting, filtering, and graph traversal FTW!

The official source code repository is at https://github.com/dib-lab/khmer and project documentation is available online at http://khmer.readthedocs.io.
See http://khmer.readthedocs.io/en/stable/introduction.html for an overview of the khmer project.

Getting help
------------

See http://khmer.readthedocs.io/en/stable/user/getting-help.html for more details, but in brief:

-  first point of contact when looking for help:
   https://github.com/dib-lab/khmer/issues
-  mailing list for **discussion**:
   http://lists.idyll.org/listinfo/khmer
-  mailing list for **announcements**:
   http://lists.idyll.org/listinfo/khmer-announce
-  email contact for project maintainers:
   khmer-project@idyll.org

Important note: cite us!
------------------------

khmer is *research software*, so you should cite us when you use it in scientific publications!
Please see the `CITATION <http://khmer.readthedocs.io/en/stable/citations.html>`__ file for citation information.

The khmer library is a project of the `Lab for Data Intensive Biology <http://ivory.idyll.org/lab/>`__ at UC Davis, and includes contributions from its members, collaborators, and friends.

Quick install
-------------

::

    pip install khmer
    pytest --pyargs khmer -m 'not known_failing and not jenkins and not huge and not linux'

See https://khmer.readthedocs.io/en/stable/user/install.html for more detailed installation instructions.

Contributing
------------

We welcome contributions to khmer from the community!
If you're interested in modifying khmer or contributing to its ongoing development see https://khmer.readthedocs.io/en/stable/dev/getting-started.html.

.. |Research software impact| image:: http://depsy.org/api/package/pypi/khmer/badge.svg
   :target: http://depsy.org/package/python/khmer
.. |Supported Python versions| image:: https://img.shields.io/pypi/pyversions/khmer.svg
.. |khmer build status| image:: https://img.shields.io/travis/dib-lab/khmer.svg
   :target: https://travis-ci.org/dib-lab/khmer
.. |Test coverage| image:: https://img.shields.io/codecov/c/github/dib-lab/khmer.svg
   :target: https://codecov.io/github/dib-lab/khmer
.. |BSD-3 licensed| image:: https://img.shields.io/badge/license-BSD%203--Clause-blue.svg
   :target: https://github.com/dib-lab/khmer/blob/master/LICENSE
