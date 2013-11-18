all:
	python setup.py build_ext -i

install:
	python setup.py install

dist:
	python setup.py sdist

clean:
	python setup.py clean --all
	cd lib && make clean
	cd tests && rm -rf khmertest_*
	rm -f khmer/_khmermodule.so

doc: FORCE
	python setup.py build_sphinx --fresh-env
	@echo ''
	@echo '--> docs in build/sphinx/html <--'
	@echo ''

cppcheck-result.xml: FORCE
	cppcheck --std=posix --platform=unix64 -j8 --enable=all -I lib/ \
		-i lib/zlib/ -i lib/bzip2/ -DVALIDATE_PARTITIONS \
		--xml lib 2> cppcheck-result.xml

cppcheck: FORCE
	cppcheck --std=posix --platform=unix64 -j8 --enable=all -I lib/ \
		-i lib/zlib/ -i lib/bzip2/ -DVALIDATE_PARTITIONS lib 

pep8: FORCE
	pep8 setup.py khmer/ scripts/ tests/

autopep8: FORCE
	autopep8 setup.py khmer/ scripts/ tests/ --recursive --in-place --pep8-passes 2000 --verbose

lib:
	cd lib && \
	$(MAKE)

test: all
	python setup.py nosetests

FORCE:
