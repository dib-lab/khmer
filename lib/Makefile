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


### NOTE: No user-serviceable parts below this line! ###

CXXFLAGS=
CXX_WARNING_FLAGS=-Wall
CXX_OPTIMIZATION_FLAGS=-O3
CXX_SHARED_LIB_FLAGS=-fPIC
CXXFLAGS+= $(CXX_WARNING_FLAGS) $(CXX_OPTIMIZATION_FLAGS) $(CXX_SHARED_LIB_FLAGS)

CFLAGS=
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

export CXX
export CXXFLAGS
export LIBS

VERSION = $(shell python get_version.py)

export VERSION

CXXFLAGS+= -DVERSION=$(VERSION)

NO_UNIQUE_RC=0
CXXFLAGS+= -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)

# Third party dependencies
ZLIB_DIR=../third-party/zlib
ZLIB_OBJS_BASE=\
       adler32.o compress.o crc32.o uncompr.o deflate.o trees.o \
       zutil.o inflate.o infback.o inftrees.o inffast.o
ZLIB_OBJS=$(addprefix $(ZLIB_DIR)/, $(ZLIB_OBJS_BASE))

BZIP2_DIR=../third-party/bzip2
BZIP2_OBJS_BASE= \
       blocksort.o huffman.o crctable.o randtable.o compress.o \
       decompress.o bzlib.o
BZIP2_OBJS=$(addprefix $(BZIP2_DIR)/, $(BZIP2_OBJS_BASE))


#DRV_PROGS=bittest # test-StreamReader test-CacheManager test-Parser test-HashTables
DRV_PROGS+=#graphtest #consume_prof
AUX_PROGS=#ht-diff


CORE_OBJS= error.o khmer_config.o thread_id_map.o trace_logger.o perf_metrics.o kmer_hash.o
PARSERS_OBJS= read_parsers.o
LIBKHMER_OBJS= $(CORE_OBJS) $(PARSERS_OBJS) hashtable.o hashbits.o subset.o counting.o scoringmatrix.o $(ZLIB_OBJS) $(BZIP2_OBJS)

KHMER_HEADERS= counting.hh error.hh hashbits.hh hashtable.hh khmer_config.hh khmer_exception.hh khmer.hh \
               kmer_hash.hh labelhash.hh perf_metrics.hh primes.hh read_aligner.hh read_parsers.hh \
	       scoringmatrix.hh subset.hh thread_id_map.hh trace_logger.hh

all: $(ZLIB_OBJS) $(BZIP2_OBJS) libkhmer test

clean:
	(cd $(ZLIB_DIR) && make distclean)
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so clean)
	rm -f *.o *.a *.so $(DRV_PROGS) $(AUX_PROGS) zlib/Makefile

test: $(DRV_PROGS) $(AUX_PROGS)

DRV_TEST_STREAM_READER_OBJS=test-StreamReader.o read_parsers.o $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_CACHE_MANAGER_OBJS=test-CacheManager.o read_parsers.o $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_PARSER_OBJS=test-Parser.o read_parsers.o $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_HASHTABLES_OBJS= \
	test-HashTables.o counting.o hashbits.o hashtable.o subset.o \
	$(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
HT_DIFF_OBJS=ht-diff.o counting.o hashtable.o $(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_COLORS_OBJS= \
	counting.o hashbits.o hashtable.o subset.o \
	$(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
test-StreamReader: $(DRV_TEST_STREAM_READER_OBJS)
	$(CXX) -o $@ $(DRV_TEST_STREAM_READER_OBJS) $(LIBS)

test-CacheManager: $(DRV_TEST_CACHE_MANAGER_OBJS)
	$(CXX) -o $@ $(DRV_TEST_CACHE_MANAGER_OBJS) $(LIBS) -fopenmp

test-Parser: $(DRV_TEST_PARSER_OBJS)
	$(CXX) -o $@ $(DRV_TEST_PARSER_OBJS) $(LIBS) -fopenmp

test-HashTables: $(DRV_TEST_HASHTABLES_OBJS)
	$(CXX) -o $@ $(DRV_TEST_HASHTABLES_OBJS) $(LIBS) -fopenmp


ht-diff: $(HT_DIFF_OBJS)
	$(CXX) -o $@ $(HT_DIFF_OBJS) $(LIBS)

bittest: bittest.o kmer_hash.o
	$(CXX) $(CXXFLAGS) -o $@ bittest.o kmer_hash.o

# NOTE: Disabled due to broken constructor call.
#graphtest: graphtest.o kmer_hash.o hashtable.o
#	$(CXX) -o $@ graphtest.o kmer_hash.o hashtable.o

# NOTE: Disabled due to broken constructor call.
#consume_prof: consume_prof.o hashtable.o kmer_hash.o $(PARSERS_OBJS)
#	$(CXX) -o $@ consume_prof.o hashtable.o kmer_hash.o $(PARSERS_OBJS) $(LIBS)

$(ZLIB_OBJS):
	(cd $(ZLIB_DIR) && ./configure && make CFLAGS+='-fPIC')

$(BZIP2_OBJS):
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so CFLAGS+='-fPIC')

%.o: %.cc $(KHMER_HEADERS)

libkhmer: libkhmer.so libkhmer.a

libkhmer.so: $(LIBKHMER_OBJS)
	$(CXX) $(CXXFLAGS) -shared -o $@ $(LIBKHMER_OBJS)

libkhmer.a: $(LIBKHMER_OBJS)
	ar cq $@ $(LIBKHMER_OBJS)
	ranlib $@

hllcounter.o: hllcounter.cc hllcounter.hh khmer.hh

test-StreamReader.o: test-StreamReader.cc read_parsers.hh

test-CacheManager.o: test-CacheManager.cc read_parsers.hh
	$(CXX) $(CXXFLAGS) -c -o $@ test-CacheManager.cc -fopenmp

test-Parser.o: test-Parser.cc read_parsers.hh
	$(CXX) $(CXXFLAGS) -c -o $@ test-Parser.cc -fopenmp

test-HashTables.o: test-HashTables.cc read_parsers.hh primes.hh
	$(CXX) $(CXXFLAGS) -c -o $@ test-HashTables.cc -fopenmp

ht-diff.o: counting.hh hashtable.hh kmer_hash.hh khmer.hh
