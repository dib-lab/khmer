all: setup-dynamic.py
	python setup.py build_ext -i

setup-dynamic.py: setup-dynamic.py.in
	sed \
		-e 's/@DEFINE_KHMER_EXTRA_SANITY_CHECKS@/$(DEFINE_KHMER_EXTRA_SANITY_CHECKS)/g' \
		-e 's/@CXX_DEBUG_FLAGS@/$(CXX_DEBUG_FLAGS)/g' \
		setup-dynamic.py.in > setup-dynamic.py

test:	all
	nosetests -v -x

clean:
	-python setup.py clean
	/bin/rm -f setup-dynamic.py
	/bin/rm -f khmer/*.so
