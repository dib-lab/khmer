#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
from nose.plugins.attrib import attr
import versioneer
import os


@attr('jenkins')
def test_python_and_c_match():
    # checks c++ compiler option version against versioneer version
    # (respectively)
    print 'c++ version {}:'.format(khmer.__version_cpp__())
    print 'versioneer (python) version: {}'.format(khmer.__version__)
    assert khmer.__version_cpp__() == khmer.__version__
