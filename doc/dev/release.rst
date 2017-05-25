..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2013-2015 Michigan State University
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

================================
Releasing a new version of khmer
================================

This document is for khmer release managers, and details the process
for making a new release of the khmer project.

How to make a khmer release candidate
-------------------------------------

Michael R. Crusoe, Luiz Irber, and C. Titus Brown have all been
release makers, following this checklist by MRC.

#. Announce a few days ahead of time that you will cut a release. This will
   slow down the rate at which PRs are being merged. When all outstanding PRs
   scheduled for this release have been merged, announce it. From now on no
   more merges to master. PRs to fix oversights or bugs follow the usual rule
   of "author can't merge". When you reach the end of this checklist announce
   it. Back to normal. If completing this checklist takes longer than a few
   hours consider allowing merges to master, and starting from the top with
   cutting a release.

#. The below should be done in a clean checkout::

        cd `mktemp -d`
        git clone git@github.com:dib-lab/khmer.git
        cd khmer

#. (Optional) Check for updates to versioneer::

        pip install --upgrade versioneer
        versioneer install
        git diff --staged

        git commit -m -a "new version of versioneer.py"
        # or
        git checkout -- versioneer.py khmer/_version.py khmer/__init__.py MANIFEST.in

#. Review the git logs since the last release and diffs (if needed) and ensure
   that ``CHANGELOG.md`` is up to date (presumably peer code review has ensured
   that it is)::

        git log --minimal --patch `git describe --tags --always --abbrev=0`..HEAD

#. Review the issue list for any new bugs that will not be fixed in this
   release. Add them to ``doc/user/known-issues.rst``

#. Check for new authors (``git log --format='%aN' v2.0... | sort -u`` lists all
   committers since the v2.0 tag). Update ``.mailmap`` to normalize their email address
   and name spelling. If they want to opt out update the ``list-*`` Makefile
   targets to exclude them. Run ``make list-citation`` and adapt the output to
   the relevant parts of ``CITATION``, ``setup.py``, ``doc/index.rst``.

#. Verify that the build is clean: https://api.travis-ci.org/dib-lab/khmer.svg?branch=master

#. Set your new version number and release candidate::

        new_version=2.2
        rc=rc1

   and then tag the release candidate with the new version number prefixed by
   the letter 'v'::

        git tag v${new_version}-${rc}
        git push --tags git@github.com:dib-lab/khmer.git

#. Test the release candidate. Bonus: repeat on Mac OS X::

        cd ..
        virtualenv testenv1
        virtualenv testenv2
        virtualenv testenv3
        virtualenv testenv4


        # First we test the tag
        cd testenv1
        source bin/activate
        git clone --depth 1 --branch v${new_version}-${rc} https://github.com/dib-lab/khmer.git
        cd khmer
        make install-dependencies
        make test
        normalize-by-median.py --version 2>&1 \
            | grep khmer\ ${new_version}-${rc} \
            && echo 1st manual version check passed
        pip uninstall -y khmer; pip uninstall -y khmer; make install
        mkdir ../not-khmer # make sure py.test executes tests
                           # from the installed khmer module
        # you might want to add 'and not huge' to the test selection
        pushd ../not-khmer; pytest --pyargs khmer.tests -m 'not known_failing'; popd


        # Secondly we test via pip
        cd ../../testenv2
        source bin/activate
        pip install -U setuptools==3.4.1
        pip install -e git+https://github.com/dib-lab/khmer.git@v${new_version}-${rc}#egg=khmer
        cd src/khmer
        make install-dependencies
        make dist
        make test
        cp dist/khmer*tar.gz ../../../testenv3/
        pip uninstall -y khmer; pip uninstall -y khmer; make install
        cd ../.. # no subdir named khmer here, safe for testing installed khmer module
        normalize-by-median.py --version 2>&1 \
            | grep khmer\ ${new_version}-${rc} \
            && echo 2nd manual version check passed
        pytest --pyargs khmer.tests -m 'not known_failing'


        # Is the distribution in testenv2 complete enough to build another
        # functional distribution?
        cd ../testenv3/
        source bin/activate
        pip install -U setuptools==3.4.1
        pip install khmer*tar.gz
        pip install pytest
        tar xzf khmer*tar.gz
        cd khmer*
        make dist
        make test
        pip uninstall -y khmer; pip uninstall -y khmer; make install
        mkdir ../not-khmer
        pushd ../not-khmer ; pytest --pyargs khmer.tests -m 'not known_failing' ; popd

#. Publish the new release on the testing PyPI server.  You will need
   to change your PyPI credentials as documented here:
   https://wiki.python.org/moin/TestPyPI.  You may need to re-register::

        python setup.py register --repository test

   Now, upload the new release::

        python setup.py sdist upload -r test

   Test the PyPI release in a new virtualenv::

        cd ../../testenv4
        source bin/activate
        pip install -U setuptools==3.4.1
        pip install screed pytest
        pip install -i https://testpypi.python.org/pypi --pre --no-clean khmer
        pytest --pyargs khmer.tests -m 'not known_failing'
        normalize-by-median.py --version 2>&1 \
            | grep khmer\ ${new_version}-${rc} \
            && echo 3rd manual version check passed
        cd build/khmer
        make test

#. Do any final acceptance tests.

#. Make sure any release notes are merged into ``doc/release-notes/``.

How to make a final release
---------------------------

When you've got a thoroughly tested release candidate, cut a release like
so:

#. Create the final tag and publish the new release on PyPI (requires an
   authorized account).::

        cd ../../../khmer
        git tag v${new_version}
        python setup.py register sdist upload

#. Delete the release candidate tag and push the tag updates to GitHub.::

        git tag -d v${new_version}-${rc}
        git push git@github.com:dib-lab/khmer.git
        git push --tags git@github.com:dib-lab/khmer.git

#. Add the release on GitHub, using the tag you just pushed.  Name
   it 'version X.Y.Z', and copy and paste in the release notes.

#. Make a binary wheel on OS X.::

        virtualenv build
        cd build
        source bin/activate
        pip install -U setuptools==3.4.1 wheel
        pip install --no-clean khmer==${new_version}
        cd build/khmer
        ./setup.py bdist_wheel upload

#. Update Read the Docs to point to the new version. Visit
   https://readthedocs.io/builds/khmer/ and 'Build Version: master' to pick up
   the new tag. Once that build has finished check the "Activate" box next to
   the new version at https://readthedocs.io/dashboard/khmer/versions/ under
   "Choose Active Versions". Finally change the default version at
   https://readthedocs.io/dashboard/khmer/advanced/ to the new version.

#. Delete any RC tags created::

        git tag -d ${new_version}-${rc}
        git push origin :refs/tags/${new_version}-${rc}

#. Tweet about the new release.

#. Send email including the release notes to khmer@lists.idyll.org
   and khmer-announce@lists.idyll.org


Setuptools Bootstrap
--------------------

`ez_setup.py` is from https://bitbucket.org/pypa/setuptools/raw/bootstrap/

Before major releases it should be examined to see if there are new
versions available and if the change would be useful


Versioning Explanation
----------------------

Versioneer, from https://github.com/warner/python-versioneer, is used to
determine the version number and is called by Setuptools and Sphinx. See the
files ``versioneer.py``, the top of ``khmer/__init__.py``,
``khmer/_version.py``, ``setup.py``, and ``doc/conf.py`` for the
implementation.

The version number is determined through several methods: see
https://github.com/warner/python-versioneer#version-identifiers

If the source tree is from a git checkout then the version number is derived by
``git describe --tags --dirty --always``. This will be in the format
``${tagVersion}-${commits_ahead}-${revision_id}-${isDirty}``. Example:
``v0.6.1-18-g8a9e430-dirty``

If from an unpacked tarball then the name of the directory is queried.

Lacking either of the two git-archive will record the version number at the top
of ``khmer/_version.py`` via the ``$Format:%d$`` and ``$Format:%H$``
placeholders enabled by the "export-subst" entry in ``.gitattributes``.

Non source distributions will have a customized ``khmer/_version.py`` that
contains hard-coded version strings. (see ``build/*/khmer/_version.py`` after a
``python setup.py build`` for an example)

``ez_setup.py`` bootstraps Setuptools (if needed) by downloading and installing
an appropriate version
