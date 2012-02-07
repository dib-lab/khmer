NO_UNIQUE_RC=0
CXXFLAGS=-g -fPIC -Wall -O3 -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)
CXXFLAGS+=-DKHMER_THREADED -fopenmp

# comment out whichever is appropriate.  can probably make this automatic ;)
SO_EXT=.so
#DYLIB_EXT=.dylib
Z_LIB_DIR=zlib-1.2.3
Z_LIB_FILES=$(Z_LIB_DIR)/*.o

all: zlib threadedParsers.o parsers.o ktable.o hashtable.o hashbits.o subset.o counting.o

clean:
	rm -f *.o $(Z_LIB_DIR)/*.o $(Z_LIB_DIR)/libz$(SO_EXT).1.2.3$(DYLIB_EXT) smpFiltering

#test: test.cc ktable.o hashtable.o

zlib: $(Z_LIB_DIR)/libz$(SO_EXT).1.2.3$(DYLIB_EXT)

$(Z_LIB_DIR)/libz$(SO_EXT).1.2.3$(DYLIB_EXT):
	cd $(Z_LIB_DIR); ./configure --shared; make; rm minigzip.o; rm example.o

smpFiltering: threadedParsers.o parsers.o counting.o hashtable.o ktable.o zlib
	$(CXX) $(CXXFLAGS) -o smpFiltering smpFiltering.cc threadedParsers.o parsers.o counting.o hashtable.o ktable.o $(Z_LIB_FILES) -fopenmp

parsetest: parsers.o 
	$(CXX) -o parsers parsers.o $(Z_LIB_FILES)

bittest: bittest.o ktable.o
	$(CXX) -o bittest bittest.o ktable.o

graphtest: graphtest.o ktable.o hashtable.o
	$(CXX) -o graphtest graphtest.o ktable.o hashtable.o

load_graph: load_graph.o hashbits.o hashtable.o threadedParsers.o parsers.o zlib
	$(CXX) $(CXXFLAGS) -pg -o load_graph load_graph.cc hashbits.o hashtable.o subset.o ktable.o threadedParsers.o parsers.o $(Z_LIB_FILES) -fopenmp

consume_prof: consume_prof.o hashtable.o ktable.o
	$(CXX) -pg -o consume_prof consume_prof.o hashtable.o ktable.o

threadedParsers.o: threadedParsers.cc threadedParsers.hh

parsers.o: parsers.cc parsers.hh

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh ktable.hh khmer.hh

intertable.o: intertable.cc intertable.hh ktable.hh khmer.hh

hashbits.o: hashbits.cc hashbits.hh subset.hh hashtable.hh ktable.hh khmer.hh counting.hh

subset.o: subset.cc subset.hh hashbits.hh ktable.hh khmer.hh

counting.o: counting.cc counting.hh hashtable.hh ktable.hh khmer.hh
