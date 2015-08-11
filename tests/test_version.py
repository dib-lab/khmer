from __future__ import print_function, unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
import khmer
from nose.plugins.attrib import attr


@attr('jenkins')
def test_python_and_c_match():
    # checks c++ compiler option version against versioneer version
    # (respectively)
    print('c++ version {0}:'.format(khmer.__version_cpp__()))
    print('versioneer (python) version: {0}'.format(khmer.__version__))
    assert khmer.__version_cpp__() == khmer.__version__


def test_python_and_c_match_base():
    # same as above but strips off the last part which can cause problems as
    # it's a hash based on git commits which can get out-of-sync too easily
    cppver = '-'.join(khmer.__version_cpp__().split('-')[0:2])
    pyver = '-'.join(khmer.__version__.split('-')[0:2])
    print('c++ version {0}'.format(cppver))
    print('python version: {0}'.format(pyver))
    print('if you are seeing this, the version compiled into your cpp')
    print('objects and your versioneer stuff is out-of-sync.')
    print('try doing: make clean; make')
    assert cppver == pyver
