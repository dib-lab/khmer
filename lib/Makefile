# Use multiple threads?
# Set this variable to true if you wish the codes to use multiple threads when they can.
WANT_THREADING=false

# Profile?
# Set this variable to true if you wish to profile the codes.
WANT_PROFILING=false

# Which profiling tool to use?
# Assuming you have TAU installed and setup properly, you can instrument codes with it to get detailed multi-threaded profiling.
# Otherwise, gprof is able to give you some information without threading info.
# Choose one of: gprof, TAU
PROFILER_OF_CHOICE=TAU

# Perform extra sanity checking?
# Set this variable to true if you wish the codes to perform extra sanity checking (to the possible detriment of performance).
WANT_EXTRA_SANITY_CHECKING=false


NO_UNIQUE_RC=0

CXXFLAGS=-g -fPIC -Wall -O3 -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)
LIBS=

ifeq ($(WANT_EXTRA_SANITY_CHECKING), true)
CXXFLAGS+=-DKHMER_EXTRA_SANITY_CHECKS
endif

ifeq ($(WANT_THREADING), true)
CXXFLAGS+=-DKHMER_THREADED -fopenmp
LIBS+=-fopenmp
endif

ifeq ($(WANT_PROFILING), true)
ifeq ($(PROFILER_OF_CHOICE), TAU)
CXX=tau_cxx.sh
endif
ifeq ($(PROFILER_OF_CHOICE), gprof)
CXXFLAGS+=-pg
endif
endif

# comment out whichever is appropriate.  can probably make this automatic ;)
SO_EXT=.so
#DYLIB_EXT=.dylib
ZLIB_DIR=zlib-1.2.3
ZLIB_OBJS=$(ZLIB_DIR)/*.o
DRV_PROGS=smpFiltering load_graph bittest ktable_test #graphtest #consume_prof

PARSER_OBJS=parsers.o
MT_PARSER_OBJS=threadedParsers.o
PARSERS_OBJS=$(PARSER_OBJS) $(MT_PARSER_OBJS)

all: zlib $(PARSERS_OBJS) ktable.o hashtable.o hashbits.o subset.o counting.o test

clean:
	(cd $(ZLIB_DIR) && make clean)
	rm -f *.o $(DRV_PROGS)

test: $(DRV_PROGS)

zlib: $(ZLIB_DIR)/libz$(SO_EXT).1.2.3$(DYLIB_EXT)

$(ZLIB_DIR)/libz$(SO_EXT).1.2.3$(DYLIB_EXT):
	cd $(ZLIB_DIR); ./configure --shared; make; rm minigzip.o; rm example.o

DRV_LOAD_GRAPH_OBJS=load_graph.o counting.o hashbits.o hashtable.o subset.o ktable.o $(PARSERS_OBJS)
DRV_SMP_FILTERING_OBJS=smpFiltering.o counting.o hashtable.o ktable.o $(PARSERS_OBJS)

smpFiltering: $(DRV_SMP_FILTERING_OBJS) zlib primes.hh
	$(CXX) -o $@ $(DRV_SMP_FILTERING_OBJS) $(ZLIB_OBJS) $(LIBS)

load_graph: $(DRV_LOAD_GRAPH_OBJS) zlib primes.hh
	$(CXX) -o $@ $(DRV_LOAD_GRAPH_OBJS) $(ZLIB_OBJS) $(LIBS)

# TODO: Move 'main' to a driver program.
#parsetest: parsers.o 
#	$(CXX) -o $@ parsers.o $(ZLIB_OBJS)

bittest: bittest.o ktable.o
	$(CXX) -o $@ bittest.o ktable.o

ktable_test: ktable_test.o ktable.o hashtable.o $(PARSERS_OBJS) zlib
	$(CXX) -o $@ ktable_test.o ktable.o hashtable.o $(PARSERS_OBJS) $(ZLIB_OBJS) $(LIBS)

# NOTE: Disabled due to broken constructor call.
#graphtest: graphtest.o ktable.o hashtable.o
#	$(CXX) -o $@ graphtest.o ktable.o hashtable.o

# NOTE: Disabled due to broken constructor call.
#consume_prof: consume_prof.o hashtable.o ktable.o $(PARSERS_OBJS)
#	$(CXX) -o $@ consume_prof.o hashtable.o ktable.o $(PARSERS_OBJS) $(LIBS)

threadedParsers.o: threadedParsers.cc threadedParsers.hh

parsers.o: parsers.cc parsers.hh

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh ktable.hh khmer.hh

hashbits.o: hashbits.cc hashbits.hh subset.hh hashtable.hh ktable.hh khmer.hh counting.hh

subset.o: subset.cc subset.hh hashbits.hh ktable.hh khmer.hh

counting.o: counting.cc counting.hh hashtable.hh ktable.hh khmer.hh
