all: lib_files python_files

clean:
	cd lib && make clean
	cd python && rm -fr build khmer/*.so

lib_files:
	cd lib && make

python_files:
	cd python && python setup.py build_ext -i

test: all
	nosetests python
