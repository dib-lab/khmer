all: lib_files python_files

clean:
	cd lib && make clean
	cd python && rm -fr build khmer/*.so

doc: FORCE
	cd doc && make html

lib_files:
	cd lib && make

python_files: lib_files
	cd python && python setup.py build_ext -i

test: all
	nosetests

FORCE:
