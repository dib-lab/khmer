all:
	./setup.py build_ext --inplace

install: FORCE
	./setup.py install

dist: FORCE
	./setup.py sdist

clean: FORCE
	./setup.py clean --all
	cd lib && ${MAKE} clean
	cd tests && rm -rf khmertest_*
	rm -f khmer/_khmermodule.so

debug:
	export CFLAGS="-pg -fprofile-arcs"; python setup.py build_ext --debug \
		--inplace

doc: FORCE
	./setup.py build_sphinx --fresh-env
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
	pep8 --exclude=_version.py setup.py khmer/ scripts/ tests/

autopep8: FORCE
	autopep8 setup.py khmer/ scripts/ tests/ --recursive --in-place \
		--pep8-passes 2000 --verbose

pylint: all FORCE
	pylint -f parseable khmer/[!_]*.py khmer/__init__.py scripts/*.py tests \
		|| true

coverage.xml: FORCE
	coverage run --branch --source=scripts,khmer -m nose --with-xunit \
		--attr=\!known_failing --processes=-1
	coverage xml

coverage-gcovr.xml: FORCE
	gcovr --root=. --branches --gcov-exclude='.*zlib.*|.*bzip2.*' --xml \
		--output=coverage-gcovr.xml

nosetests.xml: all
	./setup.py nosetests --with-xunit

doxygen: FORCE
	mkdir -p doc/doxygen
	doxygen

lib:
	cd lib && \
	$(MAKE)

test: all
	./setup.py nosetests

FORCE:
