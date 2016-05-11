all: count-stuff

clean:
	-rm count-stuff 

count-stuff: count-stuff.cc
	$(CXX) -I ../../lib count-stuff.cc -o count-stuff ../../lib/liboxli.a
