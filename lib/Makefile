CXXFLAGS=-g -fPIC -O2

all: test

clean:
	rm -f *.o

test: test.cc seqfuncs.o ktable.o hashtable.o

seqfuncs.o: seqfuncs.cc seqfuncs.hh

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh
