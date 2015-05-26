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

# Compile with performance metrics turned on?
# Set this variable to true if you want to use instrumentation provided
# in the sources for performance measurement purposes
# and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_METRICS=false


PREFIX=/usr/local

### NOTE: No user-serviceable parts below this line! ###

INCLUDES= -I ../third-party/seqan/core/include/ 	\
	  -I ../third-party/zlib/ 			\
	  -I ../third-party/bzip2/ 			\
	  -I ../third-party/smhasher/

CXXFLAGS=$(INCLUDES)
CXX_WARNING_FLAGS=-Wall
CXX_OPTIMIZATION_FLAGS=-O3
CXX_SHARED_LIB_FLAGS=-fPIC
CXXFLAGS+= 				\
	   $(CXX_WARNING_FLAGS) 	\
	   $(CXX_OPTIMIZATION_FLAGS)	\
	   $(CXX_SHARED_LIB_FLAGS)

CFLAGS=$(INCLUDES)
C_WARNING_FLAGS=-Wall
C_OPTIMIZATION_FLAGS=-O3
C_SHARED_LIB_FLAGS=-fPIC
CFLAGS+= $(C_WARNING_FLAGS) $(C_OPTIMIZATION_FLAGS) $(C_SHARED_LIB_FLAGS)

LIBS=

ifeq ($(WANT_DEBUGGING), true)
CXX_DEBUG_FLAGS=-g
CXXFLAGS+= $(CXX_DEBUG_FLAGS)
CFLAGS+= $(CXX_DEBUG_FLAGS)
else
CXX_DEBUG_FLAGS=
endif

ifeq ($(WANT_EXTRA_SANITY_CHECKING), true)
DEFINE_KHMER_EXTRA_SANITY_CHECKS=-DKHMER_EXTRA_SANITY_CHECKS
CXXFLAGS+= $(DEFINE_KHMER_EXTRA_SANITY_CHECKS)
CFLAGS+= $(DEFINE_KHMER_EXTRA_SANITY_CHECKS)
else
DEFINE_KHMER_EXTRA_SANITY_CHECKS=
endif

ifeq ($(WANT_PROFILING), true)
ifeq ($(PROFILER_OF_CHOICE), TAU)
CXX=tau_cxx.sh
endif
ifeq ($(PROFILER_OF_CHOICE), gprof)
PROFILING_LIBS=-pg
CXXFLAGS+= -pg
LIBS+= $(PROFILING_LIBS)
endif
endif

ifeq ($(WITH_INTERNAL_METRICS), true)
CXXFLAGS+= -DWITH_INTERNAL_METRICS
endif

# Place POSIX threads last in linking order, if needed.
ifneq ($(shell uname), Linux)
LIBS+= -pthread
endif


HAVE_OPENMP=$(shell 						\
	      $(CXX) -fopenmp -o chkomp .check_openmp.cc 	\
	      2>/dev/null && echo true || echo false; 		\
	      rm -f chkomp)

ifeq ($(HAVE_OPENMP), true)
	CFLAGS += -fopenmp
	CXXFLAGS += -fopenmp
endif

VERSION = $(shell python get_version.py)

LIBVERSION = $(shell python get_version.py | sed -e 's/^\([^-]*\)-.*/\1/')
LIBKHMERSO=libkhmer.so.$(LIBVERSION)

CXXFLAGS+= -DVERSION=$(VERSION)

NO_UNIQUE_RC=0
CXXFLAGS+= -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)

export CXX
export CFLAGS
export CXXFLAGS
export LIBS
export VERSION


#### Third party dependencies ####
# ZLIB
ZLIB_DIR=../third-party/zlib
ZLIB_OBJS_BASE=\
	adler32.o \
	crc32.o \
	deflate.o \
	infback.o \
	inffast.o \
	inflate.o \
	inftrees.o \
	trees.o \
	zutil.o \
	compress.o \
	uncompr.o \
	gzclose.o \
	gzlib.o \
	gzread.o \
	gzwrite.o
ZLIB_PIC_OBJS_BASE=\
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
ZLIB_PIC_OBJS=$(addprefix $(ZLIB_DIR)/, $(ZLIB_PIC_OBJS_BASE))

# BZ2
BZIP2_DIR=../third-party/bzip2
BZIP2_OBJS_BASE= \
	blocksort.o \
	huffman.o \
	crctable.o \
	randtable.o \
	compress.o \
	decompress.o \
	bzlib.o

BZIP2_OBJS=$(addprefix $(BZIP2_DIR)/, $(BZIP2_OBJS_BASE))


#### khmer proper below here ####

LIBKHMER_OBJS= \
	counting.o \
	hashbits.o \
	hashtable.o \
	hllcounter.o \
	kmer_hash.o \
	labelhash.o \
	perf_metrics.o \
	read_aligner.o \
	read_parsers.o \
	subset.o \
	trace_logger.o \
	murmur3.o \
	$(BZIP2_OBJS) \
	$(ZLIB_PIC_OBJS)

KHMER_HEADERS= \
	counting.hh \
	hashbits.hh \
	hashtable.hh \
	khmer_exception.hh \
	khmer.hh \
	kmer_hash.hh \
	labelhash.hh \
	perf_metrics.hh \
	primes.hh \
	read_aligner.hh \
	read_parsers.hh \
	subset.hh \
	trace_logger.hh

TEST_PROGS = test-Colors test-read-aligner test-compile

# START OF RULES #

# The all rule comes first!
all: $(LIBKHMERSO) libkhmer.a khmer.pc

clean:
	rm -f *.o *.a *.so* khmer.pc $(TEST_PROGS)
	(cd $(ZLIB_DIR) && make distclean)
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so clean)

test: $(TEST_PROGS)

install: $(LIBKHMERSO) libkhmer.a khmer.pc $(KHMER_HEADERS)
	mkdir -p $(PREFIX)/lib $(PREFIX)/lib/pkgconfig $(PREFIX)/include/
	cp -r $(KHMER_HEADERS) \
		../third-party/smhasher/MurmurHash3.h \
		$(PREFIX)/include/
	cp khmer.pc $(PREFIX)/lib/pkgconfig/
	cp $(LIBKHMERSO) libkhmer.a $(PREFIX)/lib
	ln -sf $(PREFIX)/lib/$(LIBKHMERSO) $(PREFIX)/lib/libkhmer.so

khmer.pc: khmer.pc.in
	sed -e 's,@prefix@,$(PREFIX),'  -e 's,@VERSION@,$(VERSION),' $< >$@

$(ZLIB_OBJS):
	(cd $(ZLIB_DIR) && ./configure && make $(ZLIB_OBJS_BASE))

$(ZLIB_PIC_OBJS):
	(cd $(ZLIB_DIR) && ./configure && make $(ZLIB_PIC_OBJS_BASE))

$(BZIP2_OBJS):
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so $(BZIP2_OBJS_BASE))

# MurMur3
murmur3.o: ../third-party/smhasher/MurmurHash3.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cc $(ZLIB_OBJS) $(ZLIB_PIC_OBJS) $(BZIP2_OBJS) $(KHMER_HEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(LIBKHMERSO): $(LIBKHMER_OBJS)
	$(CXX) $(CXXFLAGS) -shared -o $@ $(LIBKHMER_OBJS)
	ln -sf $(LIBKHMERSO) libkhmer.so

libkhmer.a: $(LIBKHMER_OBJS)
	ar rcs $@ $(LIBKHMER_OBJS)
	ranlib $@

# catch-all rule for test drivers
test-%: test-%.cc libkhmer.a
	$(CXX) $(CXXFLAGS) -I . -o $@ $< libkhmer.a
