
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

from distribute_setup import use_setuptools
use_setuptools()

from distutils.core import setup
from distutils.extension import Extension

from os.path import (
    pardir	    as path_pardir,
    join	    as path_join,
)

# Execute dynamically-generated portion of setup script in current namespace.
# Please see the dynamically-generated file to learn which variables are added
# to the current namespace.
execfile( "setup-dynamic.py" )

extra_objs = [ ]
extra_objs.extend( map(
    lambda bn: path_join( path_pardir, "lib", bn + ".o" ),
    [ 
	"khmer_config", "thread_id_map", "trace_logger", "perf_metrics", 
	"read_parsers", 
	"ktable", "hashtable", "hashbits", "counting", "subset",
        "aligner", "scoringmatrix", "node", "kmer", "outline_index"
    ]
) )
extra_objs.extend( map(
    lambda bn: path_join( path_pardir, "lib", "zlib", bn + ".o" ),
    [ 
	"adler32", "compress", "crc32", "deflate", "gzio", 
	"infback", "inffast", "inflate", "inftrees", "trees", "uncompr", 
	"zutil",
    ]
) )
extra_objs.extend( map(
    lambda bn: path_join( path_pardir, "lib", "bzip2", bn + ".o" ),
    [
	"blocksort", "huffman", "crctable", "randtable", "compress",
	"decompress", "bzlib",
    ]
) )

build_depends = extra_objs
build_depends.extend( map(
    lambda bn: path_join( path_pardir, "lib", bn + ".hh" ),
    [
	"storage", "khmer", "khmer_config", "ktable", "hashtable", "counting",
        "outline_index",
    ]
) )

extension_mod_DICT = \
    {
	"sources": [ "_khmermodule.cc" ],
	"extra_compile_args": extra_compile_args,
	"extra_link_args": filter( None, [ ] ),
	"include_dirs": [ path_join( path_pardir, "lib" ), ],
	"library_dirs": [ path_join( path_pardir, "lib" ), ],
	"extra_objects": extra_objs,
	"depends": build_depends,
	"language": "c++",
	"libraries": [ "stdc++" ],
    }

extension_mod = Extension( "khmer._khmermodule", **extension_mod_DICT )
# python modules: only 'khmer'
py_mod = 'khmer'
setup_metadata = \
    {
	"name": "khmer", "version": "0.4",
	"description": 'khmer k-mer counting library',
	"author": 'C. Titus Brown, Jason Pell, and Adina Howe',
	"author_email": 'ctb@msu.edu',
	"url": 'http://ged.msu.edu/',
	"license": 'New BSD License',
	"packages": [ py_mod, ],
	"ext_modules": [ extension_mod, ],
    }

setup( **setup_metadata )

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
