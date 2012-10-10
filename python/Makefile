all: setup.py
	python setup.py build_ext -i

setup.py: setup.py.in
	sed \
		-e 's/@DEFINE_KHMER_THREADED@/$(DEFINE_KHMER_THREADED)/g' \
		-e 's/@DEFINE_KHMER_EXTRA_SANITY_CHECKS@/$(DEFINE_KHMER_EXTRA_SANITY_CHECKS)/g' \
		-e 's/@CXX_DEBUG_FLAGS@/$(CXX_DEBUG_FLAGS)/g' \
		-e 's/@CXX_THREADING_FLAGS@/$(CXX_THREADING_FLAGS)/g' \
		-e 's/@THREADING_LIBS@/$(THREADING_LIBS)/g' \
		-e 's/@CYTHON_ENABLED_BOOL@/$(CYTHON_ENABLED_BOOL)/g' \
		setup.py.in > setup.py

test:	all
	nosetests -v -x

clean:
	-python setup.py clean
	/bin/rm -f setup.py
	/bin/rm -f khmer/*.so
