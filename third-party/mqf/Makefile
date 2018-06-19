TARGETS=gqf.o	utils.o
TESTFILES = tests/CountingTests.o tests/HighLevelFunctionsTests.o tests/IOTests.o tests/tagTests.o

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-Ofast
endif

OS := $(shell uname)
ifeq ($(OS), Darwin)
	ARCH=
else
	#ARCH=-msse4.2 -D__SSE4_2_
	ARCH=
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

#
# CXX = g++ -std=c++11
# CC = g++ -std=c++11
# LD= g++ -std=c++11

CXXFLAGS :=  -std=c++11 -fPIC -Wall $(DEBUG) $(PROFILE) $(OPT)  -m64 -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function

LDFLAGS := $(LDFLAGS) $(DEBUG) $(PROFILE) $(OPT)


#
# declaration of dependencies
#

all: $(TARGETS)



# dependencies between programs and .o files

main:	main.o	gqf.o	utils.o
	$(LD) $^ $(LDFLAGS) -o $@
# dependencies between .o files and .h files

libgqf.so: gqf.o utils.o
	$(LD) $^ $(LDFLAGS) --shared -o $@

libgqf.a: gqf.o utils.o
	ar rcs $@ $^
	ranlib $@


test:  $(TESTFILES) gqf.cc test.o utils.o
	$(LD) $(LDFLAGS) -DTEST -o mqf_test test.o utils.o $(TESTFILES) gqf.c

main.o: 								 									gqf.h

# dependencies between .o files and .cc (or .c) files


gqf.o: gqf.cc gqf.h

#
# generic build rules
#







%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.c
	$(CC) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@




clean:
	rm -f *.o $(TARGETS) $(TESTS) $(TESTFILES) libgqf.a libgqf.so
