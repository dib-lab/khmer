CXXFLAGS=-g -fPIC

all: test

clean:
	rm -f *.o

test: test.cc ktable.o hashtable.o

ktable.o: ktable.cc ktable.hh

hashtable.o: hashtable.cc hashtable.hh
