TARGETS=gqf.o

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-Ofast
endif

# using the SSE4.2 extensions doesn't seem to speed things up by a large
# amount, so disabling it for the time being till we have a good automated
# test for detecting whether or not the CPU supports it
OS := $(shell uname)
ifeq ($(OS), Darwin)
	ARCH=-m64
else
	#ARCH=-msse4.2 -D__SSE4_2_
	ARCH=
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

ifneq ($(INCLUDE),)
	INCLUDE := -I$(INCLUDE)
endif

CXX = g++
CC = g++
LD= g++

CXXFLAGS =  -std=c++11 -lstdc++ -stdlib=libstdc++  -arch x86_64  -fPIC $(ARCH) -Wall $(DEBUG) $(PROFILE) $(OPT)  -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function


LDFLAGS = $(DEBUG) $(PROFILE) $(OPT)

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files



# dependencies between .o files and .cc (or .c) files

%.o: %.cc
gqf.o: gqf.cc gqf.hh

#
# generic build rules
#

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

clean:
	rm -f *.o $(TARGETS)
