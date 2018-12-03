#! /usr/bin/env python
# vim: set fileencoding=utf-8
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
"""Setup for khmer project."""

import glob
import os
import sys
from os import listdir as os_listdir
from os.path import join as path_join
from os.path import splitext
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import codecs

from skbuild import setup

import versioneer

CMDCLASS = versioneer.get_cmdclass()

SCRIPTS = []
SCRIPTS.extend([path_join("scripts", script)
                for script in os_listdir("scripts")
                if script.endswith(".py")])

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
if "-rc" in versioneer.get_version():
    CLASSIFIERS.append("Development Status :: 4 - Beta")
else:
    CLASSIFIERS.append("Development Status :: 5 - Production/Stable")


# This sorts the author list by first name rather than last name. Not worth
#     fixing for PyPI in my opinion. The sort-authors-list.py handles it
#     correctly for the citation information, but this requires a non-standard
#     library that we don't want to add as a dependency for `setup.py`.
#     -- Daniel Standage, 2017-05-21
with codecs.open('authors.csv', 'r', encoding="utf-8") as csvin:
    authors = csvin.readlines()
authors = [a.strip().split(',') for a in authors]
authorstr = ', '.join([row[0] for row in authors])
authorstr = 'Daniel Standage, ' + authorstr + ', C. Titus Brown'

SETUP_METADATA = \
    {
        "name": "khmer",
        "version": versioneer.get_version(),
        "description": 'khmer k-mer counting library',
        "long_description": open("README.rst").read(),
        "author": authorstr,
        "author_email": 'khmer-project@idyll.org',
        # "maintainer": 'Daniel Standage', # this overrides the author field
        # "maintainer_email": 'daniel.standage@gmail.com', # so don't include
        # http://docs.python.org/2/distutils/setupscript.html
        # additional-meta-data note #3
        "url": 'https://khmer.readthedocs.io/',
        "packages": ['khmer', 'khmer.tests', 'oxli', 'khmer._oxli'],
        "package_data": {'khmer/_oxli': ['*.pxd']},
        "package_dir": {'khmer.tests': 'tests'},
        "install_requires": ['screed >= 1.0', 'bz2file', 'Cython>=0.25.2'],
        "setup_requires": ["pytest-runner>=2.0,<3dev", "setuptools>=18.0",
                           "Cython>=0.25.2"],
        "extras_require": {':python_version=="2.6"': ['argparse>=1.2.1'],
                           'docs': ['sphinx', 'sphinxcontrib-autoprogram'],
                           'tests': ['pytest>=2.9'],
                           'read_aligner_training': ['simplesam']},
        "scripts": SCRIPTS,
        # "entry_points": { # Not ready for distribution yet.
        #    'console_scripts': [
        #        "oxli = oxli:main"
        #    ]
        # },
        # "platforms": '', # empty as is conveyed by the classifiers below
        # "license": '', # empty as is conveyed by the classifier below
        "include_package_data": True,
        "zip_safe": False,
        "classifiers": CLASSIFIERS,
        "python_requires": '>=3.4'
    }


setup(cmdclass=CMDCLASS, **SETUP_METADATA)
