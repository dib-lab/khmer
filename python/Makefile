all:
	sed \
		-e 's/@DEFINE_KHMER_THREADED@/$(DEFINE_KHMER_THREADED)/g' \
		-e 's/@DEFINE_KHMER_EXTRA_SANITY_CHECKS@/$(DEFINE_KHMER_EXTRA_SANITY_CHECKS)/g' \
		-e 's/@CXX_DEBUG_FLAGS@/$(CXX_DEBUG_FLAGS)/g' \
		-e 's/@CXX_THREADING_FLAGS@/$(CXX_THREADING_FLAGS)/g' \
		-e 's/@THREADING_LIBS@/$(THREADING_LIBS)/g' \
		setup.py.in > setup.py
	python setup.py build_ext -i

test:	all
	nosetests

clean:
	python setup.py clean
	/bin/rm -f khmer/*.so
