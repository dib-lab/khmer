NO_UNIQUE_RC=0
CXXFLAGS=-g -fPIC -Wall -O2 -DNO_UNIQUE_RC=$(NO_UNIQUE_RC)

Z_LIB_DIR=zlib-1.2.3
Z_LIB_FILES=$(Z_LIB_DIR)/*.o

all: zlib parsers.o ktable.o hashtable.o intertable.o hashbits.o subset.o \
	counting.o

clean:
	rm -f *.o $(Z_LIB_DIR)/*.o $(Z_LIB_DIR)/libz.so.1.2.3

#test: test.cc ktable.o hashtable.o

zlib: $(Z_LIB_DIR)/libz.so.1.2.3

$(Z_LIB_DIR)/libz.so.1.2.3:
	cd $(Z_LIB_DIR); ./configure --shared; make; rm minigzip.o; rm example.o

parsetest: parsers.o 
	$(CXX) -o parsers parsers.o $(Z_LIB_FILES)

bittest: bittest.o ktable.o
	$(CXX) -o bittest bittest.o ktable.o

graphtest: graphtest.o ktable.o hashtable.o
	$(CXX) -o graphtest graphtest.o ktable.o hashtable.o

consume_prof: consume_prof.o hashtable.o ktable.o
	$(CXX) -pg -o consume_prof consume_prof.o hashtable.o ktable.o

parsers.o: parsers.cc parsers.hh

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh ktable.hh khmer.hh

intertable.o: intertable.cc intertable.hh ktable.hh khmer.hh

hashbits.o: hashbits.cc hashbits.hh subset.hh hashtable.hh ktable.hh khmer.hh

subset.o: subset.cc subset.hh hashbits.hh ktable.hh khmer.hh

counting.o: counting.cc counting.hh hashtable.hh ktable.hh khmer.hh
