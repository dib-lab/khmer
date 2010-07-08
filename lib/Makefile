CXXFLAGS=-g -fPIC -O2 -Wall -pg

all: ktable.o hashtable.o

clean:
	rm -f *.o

#test: test.cc ktable.o hashtable.o

bittest: bittest.o ktable.o
	$(CXX) -o bittest bittest.o ktable.o

consume_prof: consume_prof.o hashtable.o ktable.o
	$(CXX) -pg -o consume_prof consume_prof.o hashtable.o ktable.o

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh
