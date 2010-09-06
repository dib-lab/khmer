CXXFLAGS=-g -fPIC -Wall -O2

all: parsers.o ktable.o hashtable.o

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
