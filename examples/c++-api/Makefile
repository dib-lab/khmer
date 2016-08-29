all: count-demo

clean:
	-rm -f count-demo

count-demo: count-demo.cc ../../lib/liboxli.a
	$(CXX) --std=c++11 -I ../../lib count-demo.cc -o count-demo ../../lib/liboxli.a

../../lib/liboxli.a:
	cd ../../lib && make
