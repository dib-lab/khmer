# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
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

# Should we use the standard system zlib and libbz2?

USE_SYSTEM_ZLIB ?= false
USE_SYSTEM_LIBBZ2 ?= false

# Profile?
# Set this variable to true if you wish to profile the codes.
WANT_PROFILING=false

# Which profiling tool to use?
# Assuming you have TAU installed and setup properly,
# you can instrument codes with it to get detailed multi-threaded profiling.
# Otherwise, gprof is able to give you some information without threading info.
# Choose one of: gprof, TAU
PROFILER_OF_CHOICE=gprof

# Perform extra sanity checking?
# Set this variable to true
# if you wish the codes to perform extra sanity checking
# (to the possible detriment of performance).
WANT_EXTRA_SANITY_CHECKING=false

# Compile with debugging symbols?
# Set this variable to true
# if you wish the codes to be built with debugging symbols
# (increases code size
# and does not always produce accurate stepping in a debugger
# when optimization is turned on).
WANT_DEBUGGING=false

PREFIX=/usr/local

### NOTE: No user-serviceable parts below this line! ###

INCLUDES= -I ../../include/ -I ../../third-party/seqan/core/include/ \
	      -I ../../third-party/smhasher/ \
				-I ../../third-party/cqf/ \
				-I ../../third-party/rollinghash

ifeq ($(USE_SYSTEM_ZLIB), false)
INCLUDES += -I ../../third-party/zlib/
endif

ifeq ($(USE_SYSTEM_LIBBZ2), false)
INCLUDES += -I ../../third-party/bzip2/
endif

# Warnings in common to C and C++
WARNINGS=-Wall

# Flags in common to C and C++
COMMON_FLAGS=-O3 -fPIC -fno-omit-frame-pointer
SEQAN_FLAGS=-DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1

# Base C/CXXFLAGS
CPPFLAGS ?=
CPPFLAGS += $(SEQAN_FLAGS)

CXXFLAGS ?=
CXXFLAGS += $(COMMON_FLAGS) $(WARNINGS)
CXXFLAGS += -Wstrict-null-sentinel -std=c++11
CXXFLAGS += $(INCLUDES) $(CPPFLAGS)

CFLAGS	 ?=
CXXFLAGS += $(COMMON_FLAGS) $(WARNINGS)
CFLAGS   += -Wshadow -Wcast-align -Wstrict-prototypes
CFLAGS   += $(INCLUDES) $(CPPFLAGS)

LDFLAGS  ?=
ifneq ($(USE_SYSTEM_ZLIB), false)
LDFLAGS  += -lz
endif

ifneq ($(USE_SYSTEM_LIBBZ2), false)
LDFLAGS  += -lbz2
endif

ifeq ($(WANT_DEBUGGING), true)
DEBUG_FLAGS=-g
CXXFLAGS += $(DEBUG_FLAGS)
CFLAGS   += $(DEBUG_FLAGS)
endif

ifeq ($(WANT_EXTRA_SANITY_CHECKING), true)
DEFINE_OXLI_EXTRA_SANITY_CHECKS=-DOXLI_EXTRA_SANITY_CHECKS
CXXFLAGS += $(DEFINE_OXLI_EXTRA_SANITY_CHECKS)
CFLAGS   += $(DEFINE_OXLI_EXTRA_SANITY_CHECKS)
endif

ifeq ($(WANT_PROFILING), true)
ifeq ($(PROFILER_OF_CHOICE), TAU)
CXX=tau_cxx.sh
endif
ifeq ($(PROFILER_OF_CHOICE), gprof)
CXXFLAGS += -pg
CFLAGS   += -pg
LDFLAGS  += -pg
endif
endif

# Place POSIX threads last in linking order, if needed.
ifneq ($(shell uname), Linux)
LDFLAGS  += -pthread
endif


HAVE_OPENMP=$(shell 						\
	      $(CXX) -fopenmp -o chkomp .check_openmp.cc 	\
	      2>/dev/null && echo true || echo false; 		\
	      rm -f chkomp)

ifeq ($(HAVE_OPENMP), true)
CXXFLAGS +=-fopenmp
CFLAGS   +=-fopenmp
endif

ifneq ($(PACKAGE_VERSION),)
VERSION = $(PACKAGE_VERSION)
else
VERSION = $(./setup.py version | grep Version | awk '{print $2}')
endif

# The ABI version of liboxli
LIB_VERSION = 1

ifeq ($(shell uname), Darwin)
SHARED_EXT   = dylib
SONAME       = liboxli.$(SHARED_EXT).$(LIB_VERSION)
SONAME_FLAGS = -install_name $(PREFIX)/lib/$(SONAME) \
	       -compatibility_version $(LIB_VERSION) \
	       -current_version $(LIB_VERSION)
else
SHARED_EXT   = so
SONAME       = liboxli.$(SHARED_EXT).$(LIB_VERSION)
SONAME_FLAGS = -Wl,-soname=$(SONAME)
endif

LIBOXLISO=$(SONAME)

CXXFLAGS += -DVERSION=$(VERSION)

NO_UNIQUE_RC=0
CXXFLAGS += -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)
CFLAGS   += -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)

export CXX
export CFLAGS
export CXXFLAGS
export LDFLAGS
export VERSION


#### Third party dependencies ####
# ZLIB, use .lo not .o, so we get -fPIC and other library-related flags
ZLIB_DIR=../../third-party/zlib
ZLIB_OBJS_BASE=\
	adler32.lo \
	crc32.lo \
	deflate.lo \
	infback.lo \
	inffast.lo \
	inflate.lo \
	inftrees.lo \
	trees.lo \
	zutil.lo \
	compress.lo \
	uncompr.lo \
	gzclose.lo \
	gzlib.lo \
	gzread.lo \
	gzwrite.lo

ZLIB_OBJS=$(addprefix $(ZLIB_DIR)/, $(ZLIB_OBJS_BASE))

# BZ2
BZIP2_DIR=../../third-party/bzip2
BZIP2_OBJS_BASE= \
	blocksort.o \
	huffman.o \
	crctable.o \
	randtable.o \
	compress.o \
	decompress.o \
	bzlib.o

BZIP2_OBJS=$(addprefix $(BZIP2_DIR)/, $(BZIP2_OBJS_BASE))

# Counting bloom filter
CQF_DIR=../../third-party/cqf
CQF_OBJS_BASE= gqf.o

CQF_OBJS=$(addprefix $(CQF_DIR)/, $(CQF_OBJS_BASE))

#### oxli proper below here ####

LIBOXLI_OBJS= \
	hashtable.o \
	hashgraph.o \
	hllcounter.o \
	kmer_hash.o \
	labelhash.o \
	traversal.o \
	read_aligner.o \
	read_parsers.o \
	subset.o \
	kmer_filters.o \
	assembler.o \
	alphabets.o \
	murmur3.o \
	storage.o

PRECOMILE_OBJS ?=
PRECLEAN_TARGS ?=

ifeq ($(USE_SYSTEM_ZLIB), false)
LIBOXLI_OBJS  += $(ZLIB_OBJS)
PRECOMILE_OBJS += $(ZLIB_OBJS)
PRECLEAN_TARGS += zlibclean
endif

ifeq ($(USE_SYSTEM_LIBBZ2), false)
LIBOXLI_OBJS  += $(BZIP2_OBJS)
PRECOMILE_OBJS += $(BZIP2_OBJS)
PRECLEAN_TARGS += libbz2clean
endif

LIBOXLI_OBJS  += $(CQF_OBJS)
PRECOMILE_OBJS += $(CQF_OBJS)
PRECLEAN_TARGS += libcqfclean

HEADERS= \
	hashtable.hh \
	hashgraph.hh \
	hllcounter.hh \
	oxli_exception.hh \
	oxli.hh \
	kmer_hash.hh \
	labelhash.hh \
	traversal.hh \
	read_aligner.hh \
	read_parsers.hh \
	subset.hh \
	kmer_filters.hh \
	assembler.hh \
	alphabets.hh \
	storage.hh
OXLI_HEADERS = $(addprefix ../../include/oxli/,$(HEADERS))

# START OF RULES #

# The all rule comes first!
all: $(LIBOXLISO) liboxli.a oxli.pc

zlibclean:
	(cd $(ZLIB_DIR) && make distclean)
libbz2clean:
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so clean)
libcqfclean:
	(cd $(CQF_DIR) && make clean)

clean: $(PRECLEAN_TARGS)
	rm -f *.o *.a *.$(SHARED_EXT)* oxli.pc $(TEST_PROGS)

install: $(LIBOXLISO) liboxli.a oxli.pc $(OXLI_HEADERS)
	rm -rf $(PREFIX)/include/oxli $(PREFIX)/include/khmer
	mkdir -p $(PREFIX)/lib $(PREFIX)/lib/pkgconfig $(PREFIX)/include/oxli
	cp -r $(OXLI_HEADERS) 				\
		../../third-party/smhasher/MurmurHash3.h 	\
		$(PREFIX)/include/oxli/
	cp oxli.pc $(PREFIX)/lib/pkgconfig/
	cp $(LIBOXLISO) liboxli.a $(PREFIX)/lib/
	ln -sf $(PREFIX)/lib/$(LIBOXLISO) $(PREFIX)/lib/liboxli.$(SHARED_EXT)

oxli.pc: oxli.pc.in
	sed -e 's,@prefix@,$(PREFIX),'  -e 's,@VERSION@,$(VERSION),' $< >$@

$(ZLIB_OBJS):
	(cd $(ZLIB_DIR) && ./configure && make $(ZLIB_OBJS_BASE))

$(BZIP2_OBJS):
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so $(BZIP2_OBJS_BASE))

$(CQF_OBJS):
	(cd $(CQF_DIR) && make)

# MurMur3
murmur3.o: ../../third-party/smhasher/MurmurHash3.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cc $(PRECOMILE_OBJS) $(OXLI_HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $<

$(LIBOXLISO): $(LIBOXLI_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SONAME_FLAGS) -shared -o $@ $^
	ln -sf $(SONAME) liboxli.$(SHARED_EXT)

liboxli.a: $(LIBOXLI_OBJS)
	ar rcs $@ $^
	ranlib $@
