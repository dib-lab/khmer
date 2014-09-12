#!/usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
""" Extracts the version of the khmer project. """

import sys
import pkg_resources

try:
    print pkg_resources.get_distribution(  # pylint: disable=E1103
        'khmer').version
except pkg_resources.DistributionNotFound:
    print 'To build the khmer library, the distribution information'
    print 'has to be available.  Either install the package into your'
    print 'development environment or run "setup.py develop" to setup the'
    print 'metadata.  A virtualenv is recommended!'
    sys.exit(1)
del pkg_resources
