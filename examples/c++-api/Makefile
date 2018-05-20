CXX=c++
CXXFLAGS=--std=c++11 \
       -I ../../include/ \
	   -I ../../third-party/smhasher \
	   -I ../../third-party/mqf \
	   -I ../../third-party/seqan/core/include/ \
	   -I ../../third-party/rollinghash

TESTS=exact-counting bloom consume

ifneq ($(INCLUDE),)
	INCLUDE := -I$(INCLUDE)
endif

all: $(TESTS)

%: %.cc ../../src/oxli/liboxli.a
	$(CXX) -v $(CXXFLAGS) ../../src/oxli/liboxli.a  $< -o $@ -lstdc++

../../src/oxli/liboxli.a:
	cd ../../src/oxli && make

run: all
	./exact-counting
	./bloom
	./consume reads.fastq

clean:
	-rm -f exact-counting bloom consume
