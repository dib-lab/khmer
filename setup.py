#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup
from setuptools import find_packages
from setuptools import Extension

from os.path import (
    join	    as path_join,
)

from os import (
    listdir	    as os_listdir
)

from subprocess import call

zlibdir = 'lib/zlib'
bzip2dir = 'lib/bzip2'

extra_objs = [ ]
extra_objs.extend( map(
    lambda bn: path_join( "lib", "zlib", bn + ".o" ),
    [ 
	"adler32", "compress", "crc32", "deflate", "gzio", 
	"infback", "inffast", "inflate", "inftrees", "trees", "uncompr", 
	"zutil",
    ]
) )
extra_objs.extend( map(
    lambda bn: path_join( "lib", "bzip2", bn + ".o" ),
    [
	"blocksort", "huffman", "crctable", "randtable", "compress",
	"decompress", "bzlib",
    ]
) )

build_depends = list(extra_objs)
build_depends.extend( map(
    lambda bn: path_join( "lib", bn + ".hh" ),
    [
	"storage", "khmer", "khmer_config", "ktable", "hashtable", "counting",
    ]
) )

sources = [ "khmer/_khmermodule.cc" ]
sources.extend( map(
	lambda bn: path_join( "lib", bn + ".cc"),
	[
		"khmer_config", "thread_id_map", "trace_logger", "perf_metrics", 
		"read_parsers", "ktable", "hashtable", "hashbits", "counting", "subset",
	        "aligner", "scoringmatrix", "node", "kmer",
	]
) )

extension_mod_DICT = \
    {
	"sources": sources,
	"extra_compile_args": [ '-Wall', '-O3', "-pg", "-fprofile-arcs", "-ftest-coverage","-O0"],
#	"extra_link_args": filter( None, [ ] ),
	"include_dirs": [ "lib", ],
	"library_dirs": [ "lib", ],
	"extra_objects": extra_objs,
	"depends": build_depends,
	"language": "c++",
	"libraries": [ "stdc++", "gcov" ],
    }

extension_mod = Extension( "khmer._khmermodule", **extension_mod_DICT )

scripts = []
scripts.extend( [ path_join("scripts", script)
	for script in os_listdir("scripts")
	if script.endswith(".py")])

setup_metadata = \
    {
	"name": "khmer",
	"version": "0.6",
	"description": 'khmer k-mer counting library',
	"long_description": open("README.md").read(),
	"author": 'Michael R. Crusoe and Greg Edvenson and Jordan Fish and Adina Howe and Eric McDonald and Joshua Nahum and Kaben Nanlohy and Jason Pell and Jared Simpson and C. S. Welcher and Qingpeng Zhang and C. Titus Brown',
	"author_email": 'ctb@msu.edu',
	"maintainer": 'Michael R. Crusoe',
        "maintainer_email": 'mcrusoe@msu.edu',
	"url": 'http://ged.msu.edu/',
	"packages": [ 'khmer' ],
	"install_requires": [ "screed >= 0.7", 'argparse >= 1.2.1', ],
	"setup_requires": [ 'nose >= 1.0', 'setuptools-git >= 0.3', ],
	"scripts": scripts,
	"ext_modules": [ extension_mod, ],
	#"platforms": 'TODO', #??
	"include_package_data": True,
	"classifiers":	[
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
from distutils.core import setup

class build_ext(_build_ext):
	"""Specialized Python extension builder."""
	
	def run(self):
		zlib_status = call('cd ' + zlibdir +' && ( test -f Makefile || ./configure --shared ) && make libz.a', shell=True)
		bzip2_status = call('cd ' + bzip2dir + ' && make -f Makefile-libbz2_so all', shell=True)
		_build_ext.run(self)



# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
