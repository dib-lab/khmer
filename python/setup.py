#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup
from setuptools import Extension

from os.path import (
    abspath         as path_abspath,
    dirname         as path_dirname,
    pardir	    as path_pardir,
    join	    as path_join,
)

from os import (
    chdir           as os_chdir,
    listdir	    as os_listdir
)

from subprocess import call

thisdir = path_abspath(path_dirname(__file__))
rootdir = path_abspath(path_join(thisdir, '..'))
zlibdir = path_join(rootdir, 'lib/zlib')
bzip2dir = path_join(rootdir, 'lib/bzip2')

# Execute dynamically-generated portion of setup script in current namespace.
# Please see the dynamically-generated file to learn which variables are added
# to the current namespace.
#execfile(path_join(thisdir, "setup-dynamic.py" ))

extra_objs = [ ]
extra_objs.extend( map(
    lambda bn: path_join( rootdir, "lib", "zlib", bn + ".o" ),
    [ 
	"adler32", "compress", "crc32", "deflate", "gzio", 
	"infback", "inffast", "inflate", "inftrees", "trees", "uncompr", 
	"zutil",
    ]
) )
extra_objs.extend( map(
    lambda bn: path_join( rootdir, "lib", "bzip2", bn + ".o" ),
    [
	"blocksort", "huffman", "crctable", "randtable", "compress",
	"decompress", "bzlib",
    ]
) )

build_depends = list(extra_objs)
build_depends.extend( map(
    lambda bn: path_join( rootdir, "lib", bn + ".hh" ),
    [
	"storage", "khmer", "khmer_config", "ktable", "hashtable", "counting",
    ]
) )

sources = [ path_join( thisdir, "_khmermodule.cc") ]
sources.extend( map(
	lambda bn: path_join( rootdir, "lib", bn + ".cc"),
	[
		"khmer_config", "thread_id_map", "trace_logger", "perf_metrics", 
		"read_parsers", "ktable", "hashtable", "hashbits", "counting", "subset",
	        "aligner", "scoringmatrix", "node", "kmer",
	]
) )

extension_mod_DICT = \
    {
	"sources": sources,
#	"extra_compile_args": extra_compile_args,
	"extra_link_args": filter( None, [ ] ),
	"include_dirs": [ path_join( rootdir, "lib" ), ],
	"library_dirs": [ path_join( rootdir, "lib" ), ],
	"extra_objects": extra_objs,
	"depends": build_depends,
	"language": "c++",
	"libraries": [ "stdc++" ],
    }

extension_mod = Extension( "khmer._khmermodule", **extension_mod_DICT )
# python modules: only 'khmer'
py_mod = 'khmer'

scripts = []
scripts.extend( [path_join(rootdir, "scripts", script)
	for script in os_listdir(path_join(rootdir, "scripts"))
	if script.endswith(".py")])

setup_metadata = \
    {
	"name": "khmer",
	"version": "0.6",
	"description": 'khmer k-mer counting library',
	#"long_description": 'TODO',
	"author": 'Michael Crusoe and Greg Edvenson and Jordan Fish and Adina Howe and Eric McDonald and Joshua Nahum and Kaben Nanlohy and Jason Pell and Jared Simpson and C. S. Welcher and Qingpeng Zhang and C. Titus Brown',
	"author_email": 'ctb@msu.edu',
        "maintainer_email": 'mcrusoe@msu.edu',
	"url": 'http://ged.msu.edu/',
	"license": 'New BSD License',
	"packages": [ py_mod, ],
	"install_requires": "screed>=0.7",
	"scripts": scripts,
	"ext_modules": [ extension_mod, ],
	#"platforms": 'TODO', #??
	"package_dir": { 'khmer': path_join(thisdir, 'khmer')}
    }

os_chdir(zlibdir)
zlib_status = call('./configure --shared', shell=True)
zlib_make_status = call('make libz.a', shell=True)

os_chdir(bzip2dir)
bzip2_status = call('make -f Makefile-libbz2_so all', shell=True)

os_chdir(thisdir)

setup( **setup_metadata )

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
