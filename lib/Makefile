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

# Compile with tracing logic turned on?
# Set this variable to true if you want to use instrumentation provided 
# in the sources for debugging purposes 
# and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_TRACING=false

# Trace state transitions?
# Set this variable to true if you want to use instrumentation which reports
# on transitions which occur between the states of the various elements of the 
# processing stack.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_STATE_CHANGES=true

# Trace busywaits?
# Set this variable to true if you want to use instrumentation which reports
# on various busywaits, such as synchronization barriers, spinlock trials, and 
# polling loops.
# Spinlock trials will only be reported if 'TRACE_SPINLOCKS' is also true.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_BUSYWAITS=true

# Trace spinlocks?
# Set this variable to true if you want to use instrumentation which reports 
# on entries into and exits from spinlocks and spinlock trials.
# Spinlock trials will only be reported if 'TRACE_BUSYWAITS' is also true.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_SPINLOCKS=true

# Trace memory copies?
# Set this variable to true if you want to use instrumentation which reports
# on the sizes of memory copies between various caches and buffers.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_MEMCOPIES=true

# Trace data?
# Set this variable to true if you want to use instrumentation which reports 
# on the pieces of data handled by various levels of the processing stack.
# WARNING! This can generate *very* large trace logs - use with caution or 
# lots of free storage.
# 'WITH_INTERNAL_TRACING' must be true for this to have effect.
TRACE_DATA=false

# Compile with performance metrics turned on?
# Set this variable to true if you want to use instrumentation provided 
# in the sources for performance measurement purposes 
# and are willing to accept the overhead such instrumentation introduces.
WITH_INTERNAL_METRICS=false


### NOTE: No user-servicable parts below this line! ###


CXXFLAGS=
CXX_WARNING_FLAGS=-Wall
CXX_OPTIMIZATION_FLAGS=-O3
CXX_SHARED_LIB_FLAGS=-fPIC
CXXFLAGS+= $(CXX_WARNING_FLAGS) $(CXX_OPTIMIZATION_FLAGS) $(CXX_SHARED_LIB_FLAGS)

LIBS=

ifeq ($(WANT_DEBUGGING), true)
CXX_DEBUG_FLAGS=-g
CXXFLAGS+= $(CXX_DEBUG_FLAGS)
else
CXX_DEBUG_FLAGS=
endif

ifeq ($(WANT_EXTRA_SANITY_CHECKING), true)
DEFINE_KHMER_EXTRA_SANITY_CHECKS=-DKHMER_EXTRA_SANITY_CHECKS
CXXFLAGS+= $(DEFINE_KHMER_EXTRA_SANITY_CHECKS)
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

ifeq ($(WITH_INTERNAL_TRACING), true)
CXXFLAGS+= -DWITH_INTERNAL_TRACING
ifeq ($(TRACE_STATE_CHANGES), true)
CXXFLAGS+= -DTRACE_STATE_CHANGES
endif
ifeq ($(TRACE_BUSYWAITS), true)
CXXFLAGS+= -DTRACE_BUSYWAITS
endif
ifeq ($(TRACE_SPINLOCKS), true)
CXXFLAGS+= -DTRACE_SPINLOCKS
endif
ifeq ($(TRACE_MEMCOPIES), true)
CXXFLAGS+= -DTRACE_MEMCOPIES
endif
ifeq ($(TRACE_DATA), true)
CXXFLAGS+= -DTRACE_DATA
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


NO_UNIQUE_RC=0

CXXFLAGS+= -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)

ZLIB_DIR=zlib
ZLIB_OBJS_BASE=\
	adler32.o compress.o crc32.o gzio.o uncompr.o deflate.o trees.o \
	zutil.o inflate.o infback.o inftrees.o inffast.o
ZLIB_OBJS=$(addprefix $(ZLIB_DIR)/, $(ZLIB_OBJS_BASE))

BZIP2_DIR=bzip2
BZIP2_OBJS_BASE= \
	blocksort.o huffman.o crctable.o randtable.o compress.o \
	decompress.o bzlib.o
BZIP2_OBJS=$(addprefix $(BZIP2_DIR)/, $(BZIP2_OBJS_BASE))

#DRV_PROGS=bittest ktable_test test-Colors # test-StreamReader test-CacheManager test-Parser test-HashTables
DRV_PROGS+=#graphtest #consume_prof
DRV_PROGS=test-Colors
AUX_PROGS=#ht-diff

CORE_OBJS= error.o khmer_config.o thread_id_map.o trace_logger.o perf_metrics.o ktable.o
PARSERS_OBJS= read_parsers.o

all: $(ZLIB_OBJS) $(BZIP2_OBJS) $(CORE_OBJS) $(PARSERS_OBJS) hashtable.o hashbits.o subset.o counting.o test aligner.o scoringmatrix.o node.o kmer.o

clean:
	-(cd $(ZLIB_DIR) && make clean)
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so clean)
	rm -f *.o $(DRV_PROGS) $(AUX_PROGS)

test: $(DRV_PROGS) $(AUX_PROGS)

DRV_TEST_STREAM_READER_OBJS=test-StreamReader.o read_parsers.o $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_CACHE_MANAGER_OBJS=test-CacheManager.o read_parsers.o $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_PARSER_OBJS=test-Parser.o read_parsers.o $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_HASHTABLES_OBJS= \
	test-HashTables.o counting.o hashbits.o hashtable.o subset.o \
	$(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
HT_DIFF_OBJS=ht-diff.o counting.o hashtable.o $(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
DRV_TEST_COLORS_OBJS= \
	test-Colors.o counting.o hashbits.o hashtable.o subset.o \
	$(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
test-StreamReader: $(DRV_TEST_STREAM_READER_OBJS)
	$(CXX) -o $@ $(DRV_TEST_STREAM_READER_OBJS) $(LIBS)

test-CacheManager: $(DRV_TEST_CACHE_MANAGER_OBJS)
	$(CXX) -o $@ $(DRV_TEST_CACHE_MANAGER_OBJS) $(LIBS) -fopenmp

test-Parser: $(DRV_TEST_PARSER_OBJS)
	$(CXX) -o $@ $(DRV_TEST_PARSER_OBJS) $(LIBS) -fopenmp

test-HashTables: $(DRV_TEST_HASHTABLES_OBJS)
	$(CXX) -o $@ $(DRV_TEST_HASHTABLES_OBJS) $(LIBS) -fopenmp
test-Colors: $(DRV_TEST_COLORS_OBJS)
	$(CXX) -o $@ $(DRV_TEST_COLORS_OBJS) $(LIBS) -fopenmp


ht-diff: $(HT_DIFF_OBJS)
	$(CXX) -o $@ $(HT_DIFF_OBJS) $(LIBS)

bittest: bittest.o ktable.o
	$(CXX) $(CXXFLAGS) -o $@ bittest.o ktable.o

#ktable_test: ktable_test.o hashtable.o subset.o $(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
#	$(CXX) -o $@ ktable_test.o hashtable.o subset.o $(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS) $(LIBS)
color_test: test-Colors.o hashtable.o subset.o $(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS)
	$(CXX) -o $@ ktable_test.o hashtable.o subset.o $(PARSERS_OBJS) $(CORE_OBJS) $(ZLIB_OBJS) $(BZIP2_OBJS) $(LIBS)

# NOTE: Disabled due to broken constructor call.
#graphtest: graphtest.o ktable.o hashtable.o
#	$(CXX) -o $@ graphtest.o ktable.o hashtable.o

# NOTE: Disabled due to broken constructor call.
#consume_prof: consume_prof.o hashtable.o ktable.o $(PARSERS_OBJS)
#	$(CXX) -o $@ consume_prof.o hashtable.o ktable.o $(PARSERS_OBJS) $(LIBS)

$(ZLIB_OBJS):
	(cd $(ZLIB_DIR) && ./configure --shared && make libz.a)

$(BZIP2_OBJS):
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so all)

error.o: error.cc error.hh

khmer_config.o: khmer_config.cc khmer_config.hh

thread_id_map.o: thread_id_map.cc thread_id_map.hh

trace_logger.o: trace_logger.cc trace_logger.hh

perf_metrics.o: perf_metrics.cc perf_metrics.hh

read_parsers.o: read_parsers.cc read_parsers.hh

aligntest: aligntest.o zlib parsers.o ktable.o hashtable.o hashbits.o subset.o counting.o aligner.o scoringmatrix.o node.o kmer.o
	$(CXX) -pg -o aligntest aligntest.o aligner.o scoringmatrix.o node.o kmer.o $(Z_LIB_FILES) parsers.o ktable.o hashtable.o hashbits.o subset.o counting.o

bugtest: bugtest.o zlib parsers.o ktable.o hashtable.o hashbits.o subset.o counting.o aligner.o scoringmatrix.o node.o kmer.o
	$(CXX) -pg -o bugtest bugtest.o aligner.o scoringmatrix.o node.o kmer.o $(Z_LIB_FILES) parsers.o ktable.o hashtable.o hashbits.o subset.o counting.o

kmer.o: kmer.cc kmer.hh

node.o: node.cc node.hh

aligner.o: aligner.cc aligner.hh

scoringmatrix.o: scoringmatrix.cc scoringmatrix.hh

parsers.o: parsers.cc parsers.hh

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh ktable.hh khmer.hh

hashbits.o: hashbits.cc hashbits.hh subset.hh hashtable.hh ktable.hh khmer.hh counting.hh

subset.o: subset.cc subset.hh hashbits.hh ktable.hh khmer.hh

counting.o: counting.cc counting.hh hashtable.hh ktable.hh khmer.hh

test-StreamReader.o: test-StreamReader.cc read_parsers.hh

test-CacheManager.o: test-CacheManager.cc read_parsers.hh
	$(CXX) $(CXXFLAGS) -c -o $@ test-CacheManager.cc -fopenmp

test-Parser.o: test-Parser.cc read_parsers.hh
	$(CXX) $(CXXFLAGS) -c -o $@ test-Parser.cc -fopenmp

test-HashTables.o: test-HashTables.cc read_parsers.hh primes.hh
	$(CXX) $(CXXFLAGS) -c -o $@ test-HashTables.cc -fopenmp

ht-diff.o: counting.hh hashtable.hh ktable.hh khmer.hh

