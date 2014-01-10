#!/usr/bin/env python
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup
from setuptools import Extension

import versioneer
versioneer.versionfile_source = 'khmer/_version.py'
versioneer.versionfile_build = 'khmer/_version.py'
versioneer.tag_prefix = 'v'  # tags are like v1.2.0
versioneer.parentdir_prefix = '.'

from os.path import (
    join as path_join,
)

from os import (
    listdir as os_listdir
)

from subprocess import call

# strip out -Wstrict-prototypes; a hack suggested by
# http://stackoverflow.com/a/9740721
# proper fix coming in http://bugs.python.org/issue1222585
# numpy has a "nicer" fix:
# https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
import os
from distutils.sysconfig import get_config_vars
(opt,) = get_config_vars('OPT')
os.environ['OPT'] = " ".join(
    flag for flag in opt.split() if flag != '-Wstrict-prototypes'
)

zlibdir = 'lib/zlib'
bzip2dir = 'lib/bzip2'

extra_objs = []
extra_objs.extend(map(
    lambda bn: path_join("lib", "zlib", bn + ".o"),
    [
        "adler32", "compress", "crc32", "deflate", "gzio",
        "infback", "inffast", "inflate", "inftrees", "trees", "uncompr",
        "zutil"
    ]
))
extra_objs.extend(map(
    lambda bn: path_join("lib", "bzip2", bn + ".o"),
    [
        "blocksort", "huffman", "crctable", "randtable", "compress",
        "decompress", "bzlib",
    ]
))

build_depends = list(extra_objs)
build_depends.extend(map(
    lambda bn: path_join("lib", bn + ".hh"),
    [
        "storage", "khmer", "khmer_config", "ktable", "hashtable", "counting", "hashbits", "labelhash",
    ]
))

sources = ["khmer/_khmermodule.cc"]
sources.extend(map(
    lambda bn: path_join("lib", bn + ".cc"),
    [
        "khmer_config", "thread_id_map", "trace_logger", "perf_metrics",
        "read_parsers", "ktable", "hashtable", "hashbits", "labelhash", "counting",
        "subset", "aligner", "scoringmatrix", "node", "kmer",  
   ]
))

extension_mod_DICT = \
    {
        "sources": sources,
        "extra_compile_args": ['-O3', ],
        "include_dirs": ["lib", ],
        "library_dirs": ["lib", ],
        "extra_objects": extra_objs,
        "depends": build_depends,
        "language": "c++",
        "libraries": ["stdc++", ],
        "define_macros": [("VERSION", versioneer.get_version()), ],
    }

extension_mod = Extension("khmer._khmermodule", **extension_mod_DICT)

scripts = []
scripts.extend([path_join("scripts", script)
                for script in os_listdir("scripts")
                if script.endswith(".py")])

setup_metadata = \
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
        #"maintainer": 'Michael R. Crusoe', # this overrides the author field
        #"maintainer_email": 'mcrusoe@msu.edu', # so don't include it
        #http://docs.python.org/2/distutils/setupscript.html
        # #additiona-meta-data note #3
        "url": 'http://ged.msu.edu/',
        "packages": ['khmer'],
        "install_requires": ["screed >= 0.7.1", 'argparse >= 1.2.1', ],
        "setup_requires": ['nose >= 1.0', 'sphinx', ],
        "scripts": scripts,
        "ext_modules": [extension_mod, ],
        #"platforms": '', # empty as is conveyed by the classifiers below
        #"license": '', # empty as is conveyed by the classifier below
        "include_package_data": True,
        "classifiers":  [
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

# Only run lib setup when needed, not on every invocation
from distutils.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
        """Specialized Python extension builder."""

        def run(self):
                call('cd ' + zlibdir + ' && ( test -f Makefile || bash'
                     ' ./configure --shared ) && make libz.a',
                     shell=True)
                call('cd ' + bzip2dir + ' && make -f Makefile-libbz2_so all',
                     shell=True)
                _build_ext.run(self)

_cmdclass = versioneer.get_cmdclass()
_cmdclass.update({'build_ext': build_ext})

setup(cmdclass=_cmdclass,
      **setup_metadata)
