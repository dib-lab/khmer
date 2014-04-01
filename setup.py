#!/usr/bin/env python
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
""" Setup for khmer project. """

import ez_setup
ez_setup.use_setuptools(version="3.4.1")

import os
import sys
from os import listdir as os_listdir
from os.path import join as path_join

from setuptools import setup
from setuptools import Extension
from setuptools.command.build_ext import build_ext as _build_ext
from distutils.spawn import spawn
from distutils.sysconfig import get_config_vars
from distutils.dist import Distribution

import versioneer
versioneer.versionfile_source = 'khmer/_version.py'
versioneer.versionfile_build = 'khmer/_version.py'
versioneer.tag_prefix = 'v'  # tags are like v1.2.0
versioneer.parentdir_prefix = '.'
CMDCLASS = versioneer.get_cmdclass()

# strip out -Wstrict-prototypes; a hack suggested by
# http://stackoverflow.com/a/9740721
# proper fix coming in http://bugs.python.org/issue1222585
# numpy has a "nicer" fix:
# https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
OPT = get_config_vars('OPT')[0]
os.environ['OPT'] = " ".join(
    flag for flag in OPT.split() if flag != '-Wstrict-prototypes'
)

# We bundle tested versions of zlib & bzip2. To use the system zlib and bzip2
# change setup.cfg or use the `--libraries z,bz2` parameter which will make our
# custom build_ext command strip out the bundled versions.

ZLIBDIR = 'lib/zlib'
BZIP2DIR = 'lib/bzip2'

EXTRA_OBJS = []
EXTRA_OBJS.extend(path_join("lib", "zlib", bn + ".lo") for bn in [
    "adler32", "compress", "crc32", "deflate", "infback", "inffast", "inflate",
    "inftrees", "trees", "uncompr", "zutil"])
EXTRA_OBJS.extend(path_join("lib", "bzip2", bn + ".o") for bn in [
    "blocksort", "huffman", "crctable", "randtable", "compress", "decompress",
    "bzlib"])

BUILD_DEPENDS = list(EXTRA_OBJS)
BUILD_DEPENDS.extend(path_join("lib", bn + ".hh") for bn in [
    "khmer", "khmer_config", "kmer_hash", "hashtable", "counting",
    "hashbits", "labelhash"])

SOURCES = ["khmer/_khmermodule.cc"]
SOURCES.extend(path_join("lib", bn + ".cc") for bn in [
    "khmer_config", "thread_id_map", "trace_logger", "perf_metrics",
    "read_parsers", "kmer_hash", "hashtable", "hashbits", "labelhash",
    "counting", "subset", "aligner", "scoringmatrix", "node", "kmer"])

EXTRA_COMPILE_ARGS = ['-O3']

if sys.platform == 'darwin':
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64'])  # force 64bit only builds

EXTENSION_MOD_DICT = \
    {
        "sources": SOURCES,
        "extra_compile_args": EXTRA_COMPILE_ARGS,
        "extra_objects": EXTRA_OBJS,
        "depends": BUILD_DEPENDS,
        "language": "c++",
        "define_macros": [("VERSION", versioneer.get_version()), ],
    }

EXTENSION_MOD = Extension("khmer._khmermodule",  # pylint: disable=W0142
                          ** EXTENSION_MOD_DICT)
SCRIPTS = []
SCRIPTS.extend([path_join("scripts", script)
                for script in os_listdir("scripts")
                if script.endswith(".py")])

SETUP_METADATA = \
    {
        "name": "khmer",
        "version": versioneer.get_version(),
        "description": 'khmer k-mer counting library',
        "long_description": open("README.rst").read(),
        "author": 'Michael R. Crusoe, Greg Edvenson, Jordan Fish,'
        ' Adina Howe, Eric McDonald, Joshua Nahum, Kaben Nanlohy,'
        ' Jason Pell, Jared Simpson, Camille Scott,'
        ' Qingpeng Zhang, and C. Titus Brown',
        "author_email": 'khmer-project@idyll.org',
        # "maintainer": 'Michael R. Crusoe', # this overrides the author field
        # "maintainer_email": 'mcrusoe@msu.edu', # so don't include it
        # http://docs.python.org/2/distutils/setupscript.html
        # additiona-meta-data note #3
        "url": 'http://ged.msu.edu/',
        "packages": ['khmer'],
        "install_requires": ["screed >= 0.7.1", 'argparse >= 1.2.1', ],
        "extras_require": {'docs': ['sphinx', 'sphinxcontrib-autoprogram'],
                           'tests': ['nose >= 1.0']},
        "scripts": SCRIPTS,
        "ext_modules": [EXTENSION_MOD, ],
        # "platforms": '', # empty as is conveyed by the classifiers below
        # "license": '', # empty as is conveyed by the classifier below
        "include_package_data": True,
        "zip_safe": False,
        "classifiers": [
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Environment :: MacOS X",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: POSIX :: Linux",
            "Operating System :: MacOS :: MacOS X",
            "Programming Language :: C",
            "Programming Language :: C++",
            "Programming Language :: Python :: 2.7",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    }


class KhmerBuildExt(_build_ext):  # pylint: disable=R0904
    """Specialized Python extension builder for khmer project.

    Only run the library setup when needed, not on every invocation.

    Also strips out the bundled zlib and bzip2 libraries if
    `--libraries z,bz2` is specified or the equivalent is in setup.cfg
    """

    def run(self):
        if "z" and "bz2" not in self.libraries:
            spawn(cmd=['bash', '-c', 'cd ' + ZLIBDIR + ' && ( test configure '
                       '-nt Makefile || bash ./configure --static ) && make '
                       '-f Makefile.pic PIC'],
                  dry_run=self.dry_run)
            spawn(cmd=['bash', '-c', 'cd ' + BZIP2DIR + ' && make -f '
                       'Makefile-libbz2_so all'],
                  dry_run=self.dry_run)
        else:
            for ext in self.extensions:
                ext.extra_objects = []

        _build_ext.run(self)

CMDCLASS.update({'build_ext': KhmerBuildExt})

_DISTUTILS_REINIT = Distribution.reinitialize_command


def reinitialize_command(self, command, reinit_subcommands):
    '''
    Monkeypatch distutils.Distribution.reinitialize_command() to match behavior
    of Distribution.get_command_obj()

    This fixes issues with 'pip install -e' and './setup.py nosetests' not
    respecting the setup.cfg configuration directives for the build_ext command
    '''
    cmd_obj = _DISTUTILS_REINIT(self, command, reinit_subcommands)
    options = self.command_options.get(command)
    if options:
        self._set_command_options(  # pylint: disable=protected-access
            cmd_obj, options)
    return cmd_obj
Distribution.reinitialize_command = reinitialize_command

# pylint: disable=W0142
setup(cmdclass=CMDCLASS,
      **SETUP_METADATA)
