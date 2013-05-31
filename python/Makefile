all: setup.py
	python setup.py build_ext -i

setup.py: setup.py.in
	sed \
		-e 's/@DEFINE_KHMER_EXTRA_SANITY_CHECKS@/$(DEFINE_KHMER_EXTRA_SANITY_CHECKS)/g' \
		-e 's/@CXX_DEBUG_FLAGS@/$(CXX_DEBUG_FLAGS)/g' \
		setup.py.in > setup.py

test:	all
	nosetests -v -x

clean:
	-python setup.py clean
	/bin/rm -f setup.py
	/bin/rm -f khmer/*.so
