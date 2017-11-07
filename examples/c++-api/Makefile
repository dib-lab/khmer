CXXFLAGS=--std=c++11 \
       -I ../../include/ \
	   -I ../../third-party/smhasher \
	   -I ../../third-party/cqf \
	   -I ../../third-party/seqan/core/include/ \
	   -I ../../third-party/rollinghash

TESTS=exact-counting bloom consume

%: %.cc ../../src/oxli/liboxli.a
	$(CXX) $(CXXFLAGS) $< ../../src/oxli/liboxli.a -o $@

../../src/oxli/liboxli.a:
	cd ../../src/oxli && make

run: all
	./exact-counting
	./bloom
	./consume reads.fastq

all: $(TESTS)

clean:
	-rm -f exact-counting bloom consume
