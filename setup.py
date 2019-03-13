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

import ez_setup

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

from setuptools import setup
from setuptools import Extension
from setuptools.command.build_ext import build_ext as _build_ext
from distutils.spawn import spawn
from distutils.sysconfig import get_config_vars
from distutils.dist import Distribution
from distutils.errors import DistutilsPlatformError

import versioneer
ez_setup.use_setuptools(version="3.4.1")

CMDCLASS = versioneer.get_cmdclass()

from setuptools import Extension as CyExtension
HAS_CYTHON = False
cy_ext = 'cpp'

# strip out -Wstrict-prototypes; a hack suggested by
# http://stackoverflow.com/a/9740721
# proper fix coming in http://bugs.python.org/issue1222585
# numpy has a "nicer" fix:
# https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
OPT = get_config_vars('OPT')[0]
os.environ['OPT'] = " ".join(
    flag for flag in OPT.split() if flag != '-Wstrict-prototypes'
)

# Checking for OpenMP support. Currently clang doesn't work with OpenMP,
# so it needs to be disabled for now.
# This function comes from the yt project:
# https://bitbucket.org/yt_analysis/yt/src/f7c75759e0395861b52d16921d8ce3ad6e36f89f/yt/utilities/lib/setup.py?at=yt


def check_for_openmp():
    """Check for OpenMP support."""
    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    exit_code = 1

    if os.name == 'nt':
        return False

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CC', 'cc')

        # Attempt to compile a test script.
        # See http://openmp.org/wp/openmp-compilers/
        filename = r'test.c'
        source = open(filename, 'wt', 1)
        source.write(
            """
            #include <omp.h>
            #include <stdio.h>
            int main() {
            #pragma omp parallel
            printf("Hello from thread %d, nthreads %d",
                    omp_get_thread_num(), omp_get_num_threads());
            }
            """
        )
        with open(os.devnull, 'w') as fnull:
            exit_code = subprocess.call([compiler, '-fopenmp', filename],
                                        stdout=fnull, stderr=fnull)

        # Clean up
        source.close()
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return exit_code == 0


def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)


def build_dir():
    return path_join("build", distutils_dir_name("temp"))

# We bundle tested versions of zlib & bzip2. To use the system zlib and bzip2
# change setup.cfg or use the `--libraries z,bz2` parameter which will make our
# custom build_ext command strip out the bundled versions.


ZLIBDIR = 'third-party/zlib'
BZIP2DIR = 'third-party/bzip2'

BUILD_DEPENDS = [path_join("include", "khmer", bn + ".hh") for bn in [
    "_cpy_khmer", "_cpy_utils", "_cpy_readparsers"
]]
BUILD_DEPENDS.extend(path_join("include", "oxli", bn + ".hh") for bn in [
    "khmer", "kmer_hash", "hashtable", "labelhash", "hashgraph",
    "hllcounter", "oxli_exception", "read_aligner", "subset", "read_parsers",
    "kmer_filters", "traversal", "assembler", "alphabets", "storage"])

SOURCES = [path_join("src", "khmer", bn + ".cc") for bn in [
    "_cpy_khmer", "_cpy_utils", "_cpy_readparsers"
]]
SOURCES.extend(path_join("src", "oxli", bn + ".cc") for bn in [
    "read_parsers", "kmer_hash", "hashtable", "hashgraph",
    "labelhash", "subset", "read_aligner",
    "hllcounter", "traversal", "kmer_filters", "assembler", "alphabets",
    "storage"])

SOURCES.extend(path_join("third-party", "smhasher", bn + ".cc") for bn in [
    "MurmurHash3"])

# Don't forget to update lib/Makefile with these flags!
EXTRA_COMPILE_ARGS = ['-O3', '-std=c++11']
EXTRA_LINK_ARGS = []

if sys.platform == 'darwin':
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64', '-mmacosx-version-min=10.9'])
    EXTRA_LINK_ARGS.append('-mmacosx-version-min=10.9')

if check_for_openmp():
    EXTRA_COMPILE_ARGS.extend(['-fopenmp'])
    EXTRA_LINK_ARGS.extend(['-fopenmp'])

CP_EXTENSION_MOD_DICT = \
    {
        "sources": SOURCES,
        "extra_compile_args": EXTRA_COMPILE_ARGS,
        "extra_link_args": EXTRA_LINK_ARGS,
        "depends": BUILD_DEPENDS,
        "include_dirs": ["include", "."],
        "language": "c++",
        "define_macros": [("VERSION", versioneer.get_version()), ],
    }

EXTENSION_MODS = [Extension("khmer._khmer", ** CP_EXTENSION_MOD_DICT)]

CY_OPTS = {
    'embedsignature': True,
    'language_level': 3,
    'c_string_type': 'unicode',
    'c_string_encoding': 'utf8'
}

for cython_ext in glob.glob(os.path.join("khmer", "_oxli",
                                         "*.{0}".format(cy_ext))):

    CY_EXTENSION_MOD_DICT = \
        {
            "sources": [cython_ext, "khmer/_oxli/oxli_exception_convert.cc"],
            "extra_compile_args": EXTRA_COMPILE_ARGS,
            "extra_link_args": EXTRA_LINK_ARGS,
            "extra_objects": [path_join(build_dir(), splitext(p)[0] + '.o')
                              for p in SOURCES],
            "depends": [],
            "include_dirs": ["include", "."],
            "language": "c++",
            "define_macros": [("VERSION", versioneer.get_version()), ]
        }

    if HAS_CYTHON:
        CY_EXTENSION_MOD_DICT['cython_directives'] = CY_OPTS

    ext_name = "khmer._oxli.{0}".format(
        splitext(os.path.basename(cython_ext))[0]
    )
    EXTENSION_MODS.append(CyExtension(ext_name, ** CY_EXTENSION_MOD_DICT))

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
        "install_requires": ['screed>=1.0', 'bz2file'],
        "setup_requires": ['setuptools>=18.0'],
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
        "ext_modules": EXTENSION_MODS,
        # "platforms": '', # empty as is conveyed by the classifiers below
        # "license": '', # empty as is conveyed by the classifier below
        "include_package_data": True,
        "zip_safe": False,
        "classifiers": CLASSIFIERS,
        "python_requires": '>=3.4'
    }


class KhmerBuildExt(_build_ext):  # pylint: disable=R0904
    """Specialized Python extension builder for khmer project.

    Only run the library setup when needed, not on every invocation.

    Also strips out the bundled zlib and bzip2 libraries if
    `--libraries z,bz2` is specified or the equivalent is in setup.cfg
    """

    def run(self):
        """Run extension builder."""
        if "%x" % sys.maxsize != '7fffffffffffffff':
            raise DistutilsPlatformError("%s require 64-bit operating system" %
                                         SETUP_METADATA["packages"])

        if sys.platform == 'darwin' and 'gcov' in self.libraries:
            self.libraries.remove('gcov')

        cqfcmd = ['bash', '-c', 'cd third-party/cqf && make']
        spawn(cmd=cqfcmd, dry_run=self.dry_run)
        for ext in self.extensions:
            ext.extra_objects.append(path_join("third-party", "cqf", "gqf.o"))

        if "z" not in self.libraries:
            zcmd = ['bash', '-c', 'cd ' + ZLIBDIR + ' && ( test Makefile -nt'
                    ' configure || bash ./configure --static ) && make -f '
                    'Makefile.pic PIC']
            spawn(cmd=zcmd, dry_run=self.dry_run)
            for ext in self.extensions:
                ext.extra_objects.extend(
                    path_join("third-party", "zlib", bn + ".lo") for bn in [
                        "adler32", "compress", "crc32", "deflate", "gzclose",
                        "gzlib", "gzread", "gzwrite", "infback", "inffast",
                        "inflate", "inftrees", "trees", "uncompr", "zutil"])

        if "bz2" not in self.libraries:
            bz2cmd = ['bash', '-c', 'cd ' + BZIP2DIR + ' && make -f '
                      'Makefile-libbz2_so all']
            spawn(cmd=bz2cmd, dry_run=self.dry_run)
            for ext in self.extensions:
                ext.extra_objects.extend(
                    path_join("third-party", "bzip2", bn + ".o") for bn in [
                        "blocksort", "huffman", "crctable", "randtable",
                        "compress", "decompress", "bzlib"])
        _build_ext.run(self)


CMDCLASS.update({'build_ext': KhmerBuildExt})

_DISTUTILS_REINIT = Distribution.reinitialize_command


def reinitialize_command(self, command, reinit_subcommands):
    """Monkeypatch the original version from distutils.

    It's supposed to match the behavior of Distribution.get_command_obj()
    This fixes issues with 'pip install -e' and './setup.py test' not
    respecting the setup.cfg configuration directives for the build_ext
    command.
    """
    cmd_obj = _DISTUTILS_REINIT(self, command, reinit_subcommands)
    options = self.command_options.get(command)
    if options:
        self._set_command_options(  # pylint: disable=protected-access
            cmd_obj, options)
    return cmd_obj


Distribution.reinitialize_command = reinitialize_command


setup(cmdclass=CMDCLASS, **SETUP_METADATA)
