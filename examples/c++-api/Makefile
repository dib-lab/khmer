CFLAGS=--std=c++11 \
	   -I ../../include \
	   -I ../../third-party/smhasher \
	   -I ../../third-party/seqan/core/include/

run: all
	./exact-counting
	./bloom
	./consume reads.fastq

all: exact-counting bloom consume

clean:
	-rm -f exact-counting bloom consume

%: %.cc liboxli_a
	$(CXX) $(CFLAGS) -o $* $*.cc ../../src/oxli/liboxli.a

liboxli_a:
	cd ../../src/oxli && make
