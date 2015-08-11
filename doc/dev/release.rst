.. vim: set filetype=rst

================================
Releasing a new version of khmer
================================

This document is for khmer release managers, and details the process
for making a new release of the khmer project.

How to make a khmer release candidate
-------------------------------------

Michael R. Crusoe, Luiz Irber, and C. Titus Brown have all been
release makers, following this checklist by MRC.

#. The below should be done in a clean checkout::

        cd `mktemp -d`
        git clone git@github.com:dib-lab/khmer.git
        cd khmer

#. (Optional) Check for updates to versioneer::

        pip install --upgrade versioneer
        versioneer-installer

        git diff
        
        ./setup.py versioneer
        git diff 
        git commit -m -a "new version of versioneer.py"
        # or
        git checkout -- versioneer.py khmer/_version.py khmer/__init__.py MANIFEST.in

#. Review the git logs since the last release and diffs (if needed) and ensure
   that the ``ChangeLog`` is up to date::

        git log --minimal --patch `git describe --tags --always --abbrev=0`..HEAD

#. Review the issue list for any new bugs that will not be fixed in this
   release. Add them to ``doc/known-issues.txt``

#. Verify that the build is clean: http://ci.ged.msu.edu/job/khmer-master/

#. Submit a build to Coverity Scan if it hasn't been done
   recently. You can get the token from
   https://gitlab.msu.edu/ged-lab/ged-internal-docs/wikis/coverity-scan
   or https://scan.coverity.com/projects/621?tab=project_settings

   ::

        virtualenv coverityenv
        source coverityenv/bin/activate
        make install-dependencies
        make clean
        cov_analysis_dir=~/src/coverity/cov-analysis-linux64-7.5.0/ make coverity-build
        COVERITY_TOKEN=${COVERITY_TOKEN} make coverity-upload

#. Set your new version number and release candidate::

        new_version=1.3
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
        normalize-by-median.py --version 2>&1 | grep khmer\ ${new_version}-${rc} && \
                echo 1st manual version check passed
        pip uninstall -y khmer; pip uninstall -y khmer; make install
        mkdir ../not-khmer # if there is a subdir named 'khmer' nosetest will execute tests
        # there instead of the installed khmer module's tests
        pushd ../not-khmer; nosetests khmer --attr '!known_failing'; popd


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
        cd ../.. # no subdir named khmer here, safe for nosetesting installed khmer module
        normalize-by-median.py --version 2>&1 | grep khmer\ ${new_version}-${rc} && \
                echo 2nd manual version check passed
        nosetests khmer --attr '!known_failing'

        # Is the distribution in testenv2 complete enough to build another
        # functional distribution?
        
        cd ../testenv3/
        source bin/activate
        pip install -U setuptools==3.4.1
        pip install khmer*tar.gz
        pip install nose
        tar xzf khmer*tar.gz
        cd khmer*
        make dist
        make test
        pip uninstall -y khmer; pip uninstall -y khmer; make install
        mkdir ../not-khmer
        pushd ../not-khmer ; nosetests khmer --attr '!known_failing' ; popd

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
        pip install screed nose
        pip install -i https://testpypi.python.org/pypi --pre --no-clean khmer
        nosetests khmer --attr '!known_failing'
        normalize-by-median.py --version 2>&1 | grep khmer\ ${new_version}-${rc} && \
                echo 3rd manual version check passed
        cd build/khmer
        make test

#. Do any final testing (BaTLab and/or acceptance tests).

#. Make sure any release notes are merged into doc/release-notes/.

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
   https://readthedocs.org/builds/khmer/ and 'Build Version: master' to pick up
   the new tag. Once that build has finished check the "Activate" box next to
   the new version at https://readthedocs.org/dashboard/khmer/versions/ under
   "Choose Active Versions". Finally change the default version at
   https://readthedocs.org/dashboard/khmer/advanced/ to the new version.

#. Delete any RC tags created::

        git tag -d ${new_version}-${rc}
        git push origin :refs/tags/${new_version}-${rc}

#. Tweet about the new release.

#. Send email including the release notes to khmer@lists.idyll.org
   and khmer-announce@lists.idyll.org

BaTLab testing
--------------

The UW-Madison Build and Test Lab provides the khmer project with a free
cross-platform testing environment.

#. Connect to their head node::

        ssh mcrusoe@submit-1.batlab.org

#. Move into the khmer directory and download a release from PyPI's main server
   or the test PyPI server::

        cd khmer/
        wget https://testpypi.python.org/packages/source/k/khmer/khmer-1.0.1-rc3.tar.gz
        vim khmer-v1.0.inputs # change the 'scp_file' to point to the release
        vim khmer-v1.0.run-spec # change 'project_version' at bottom
        nmi_submit khmer-v1.0.run-spec

Setuptools Bootstrap
--------------------

ez_setup.py is from https://bitbucket.org/pypa/setuptools/raw/bootstrap/

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
