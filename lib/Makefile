CXXFLAGS=-g -fPIC -Wall -O2

all: parsers.o ktable.o hashtable.o intertable.o hashbits.o subset.o

clean:
	rm -f *.o

#test: test.cc ktable.o hashtable.o

parsers.o: parsers.cc parsers.hh

bittest: bittest.o ktable.o
	$(CXX) -o bittest bittest.o ktable.o

graphtest: graphtest.o ktable.o hashtable.o
	$(CXX) -o graphtest graphtest.o ktable.o hashtable.o

consume_prof: consume_prof.o hashtable.o ktable.o
	$(CXX) -pg -o consume_prof consume_prof.o hashtable.o ktable.o

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh


intertable.o: intertable.cc intertable.hh

hashbits.o: hashbits.cc hashbits.hh subset.hh

subset.o: subset.cc subset.hh hashbits.hh


