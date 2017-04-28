CFLAGS=--std=c++11 \
       -I ../../lib \
	   -I ../../third-party/smhasher \
	   -I ../../third-party/seqan/core/include/

run: all
	./exact-counting
	./bloom
	./consume reads.fastq

all: exact-counting bloom consume

clean:
	-rm -f exact-counting bloom consume

%: %.cc ../../lib/liboxli.a
	$(CXX) $(CFLAGS) -o $* $*.cc ../../lib/liboxli.a

../../lib/liboxli.a:
	cd ../../lib && make
