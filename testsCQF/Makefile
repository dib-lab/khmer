TARGETS= mergeTest_SameSize mergeTest_DifferentSize mergeTest_Resize mergeTest_Qbits


DEBUG=-g
OPT=
ARCH=


PROFILE=

CXX = g++ -std=c++11
CC = g++ -std=c++11
LD= g++ -std=c++11

CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I../third-party/cqf/ -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lm

#
# declaration of dependencies
#

all:	$(TARGETS)

# dependencies between programs and .o files

mergeTest_SameSize:	mergeTest_SameSize.o
mergeTest_DifferentSize: mergeTest_DifferentSize.o
mergeTest_Resize:	mergeTest_Resize.o
mergeTest_Qbits:	mergeTest_Qbits.o
# dependencies between .o files and .h files


# dependencies between .o files and .cc (or .c) files



#
# generic build rules
#

$(TARGETS):
	$(LD) $^ ../third-party/cqf/gqf.o $(LDFLAGS) -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.c
	$(CC) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

clean:
	rm -f *.o $(TARGETS)
