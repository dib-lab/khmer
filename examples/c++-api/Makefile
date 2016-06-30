all: count-stuff

clean:
	-rm -f count-stuff

count-stuff: count-stuff.cc ../../lib/liboxli.a
	$(CXX) --std=c++11 -I ../../lib count-stuff.cc -o count-stuff ../../lib/liboxli.a

../../lib/liboxli.a:
	cd ../../lib && make
