CXXFLAGS=--std=c++11 \
       -I ../../include/ \
	   -I ../../third-party/smhasher \
	   -I ../../third-party/mqf \
	   -I ../../third-party/seqan/core/include/ \
	   -I ../../third-party/rollinghash\
		 -Wl,-rpath=../../src/oxli/

TESTS=exact-counting bloom consume


ifneq ($(INCLUDE),)
	INCLUDE := -I$(INCLUDE)
endif

all: $(TESTS)

%: %.cc ../../src/oxli/liboxli.so
	$(CXX) -L../../src/oxli/liboxli.so  -lstdc++ $(CXXFLAGS) $<  -o $@ -L../..//src/oxli liboxli.so.1 -lstdc++

../../src/oxli/liboxli.so:
	cd ../../src/oxli && make

run: all
	./exact-counting
	./bloom
	./consume reads.fastq

clean:
	-rm -f exact-counting bloom consume
