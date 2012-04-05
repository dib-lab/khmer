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

DRV_PROGS=smpFiltering load_graph bittest ktable_test #graphtest #consume_prof

PARSER_OBJS=parsers.o
MT_PARSER_OBJS=threadedParsers.o
PARSERS_OBJS=$(PARSER_OBJS) $(MT_PARSER_OBJS)

all: $(ZLIB_OBJS) $(BZIP2_OBJS) $(PARSERS_OBJS) khmer_config.o ktable.o hashtable.o hashbits.o subset.o counting.o test

clean:
	(cd $(ZLIB_DIR) && make clean)
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so clean)
	rm -f *.o $(DRV_PROGS)

test: $(DRV_PROGS)

DRV_LOAD_GRAPH_OBJS=load_graph.o counting.o hashbits.o hashtable.o subset.o ktable.o khmer_config.o $(PARSERS_OBJS)
DRV_SMP_FILTERING_OBJS=smpFiltering.o counting.o hashtable.o ktable.o khmer_config.o $(PARSERS_OBJS)

smpFiltering: $(DRV_SMP_FILTERING_OBJS) $(ZLIB_OBJS) primes.hh
	$(CXX) -o $@ $(DRV_SMP_FILTERING_OBJS) $(ZLIB_OBJS) $(LIBS)

load_graph: $(DRV_LOAD_GRAPH_OBJS) $(ZLIB_OBJS) primes.hh
	$(CXX) -o $@ $(DRV_LOAD_GRAPH_OBJS) $(ZLIB_OBJS) $(LIBS)

# TODO: Move 'main' to a driver program.
#parsetest: parsers.o 
#	$(CXX) -o $@ parsers.o $(ZLIB_OBJS)

bittest: bittest.o ktable.o
	$(CXX) -o $@ bittest.o ktable.o

ktable_test: ktable_test.o ktable.o hashtable.o $(PARSERS_OBJS) $(ZLIB_OBJS)
	$(CXX) -o $@ ktable_test.o ktable.o hashtable.o $(PARSERS_OBJS) $(ZLIB_OBJS) $(LIBS)

# NOTE: Disabled due to broken constructor call.
#graphtest: graphtest.o ktable.o hashtable.o
#	$(CXX) -o $@ graphtest.o ktable.o hashtable.o

# NOTE: Disabled due to broken constructor call.
#consume_prof: consume_prof.o hashtable.o ktable.o $(PARSERS_OBJS)
#	$(CXX) -o $@ consume_prof.o hashtable.o ktable.o $(PARSERS_OBJS) $(LIBS)

$(ZLIB_OBJS):
	(cd $(ZLIB_DIR) && ./configure --shared && make libz.so.1.2.3)

$(BZIP2_OBJS):
	(cd $(BZIP2_DIR) && make -f Makefile-libbz2_so all)

khmer_config.o: khmer_config.cc khmer_config.hh

threadedParsers.o: threadedParsers.cc threadedParsers.hh

parsers.o: parsers.cc parsers.hh

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh ktable.hh khmer.hh

hashbits.o: hashbits.cc hashbits.hh subset.hh hashtable.hh ktable.hh khmer.hh counting.hh

subset.o: subset.cc subset.hh hashbits.hh ktable.hh khmer.hh

counting.o: counting.cc counting.hh hashtable.hh ktable.hh khmer.hh
