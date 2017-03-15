CFLAGS=--std=c++11 \
       -I ../../lib \
	   -I ../../third-party/smhasher \
	   -I ../../third-party/seqan/core/include/

run: all
	./exact-counting
	./bloom

all: exact-counting bloom

clean:
	-rm -f exact-counting bloom

bloom: bloom.cc ../../lib/liboxli.a
	$(CXX) $(CFLAGS) bloom.cc -o bloom ../../lib/liboxli.a

exact-counting: exact-counting.cc ../../lib/liboxli.a
	$(CXX) $(CFLAGS) exact-counting.cc -o exact-counting ../../lib/liboxli.a

../../lib/liboxli.a:
	cd ../../lib && make
