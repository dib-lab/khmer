all:
	python setup.py build_ext -i

test:	all
	nosetests

clean:
	python setup.py clean
	/bin/rm -f khmer/*.so
