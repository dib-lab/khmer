# make pep8 to check for basic Python code compliance
# make autopep8 to fix most pep8 errors
# make pylint to check Python code for enhanced compliance including naming
#  and documentation

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
	ls lib/*.cc khmer/_khmermodule.cc | grep -v test | cppcheck -DNDEBUG \
		-DVERSION=0.0.cppcheck -UNO_UNIQUE_RC --enable=all --file-list=- -j8 \
		--platform=unix64 --std=posix --xml --xml-version=2 2> cppcheck-result.xml

cppcheck: FORCE
	ls lib/*.cc khmer/_khmermodule.cc | grep -v test | cppcheck -DNDEBUG \
		-DVERSION=0.0.cppcheck -UNO_UNIQUE_RC --enable=all --file-list=- -j8 \
		--platform=unix64 --std=posix --quiet

pep8: FORCE
	pep8 --exclude=_version.py setup.py khmer/ scripts/ tests/

autopep8: FORCE
	autopep8 --recursive --in-place --exclude _version.py --ignore E309 setup.py \
		khmer/ scripts/ tests/

pylint: all FORCE
	pylint -f parseable khmer/[!_]*.py khmer/__init__.py scripts/*.py tests \
		|| true

coverage.xml: FORCE
	coverage run --branch --source=scripts,khmer -m nose --with-xunit \
		--attr=\!known_failing --processes=0
	coverage xml

coverage.html: FORCE
	coverage run --branch --source=scripts,khmer -m nose --with-xunit \
		--attr=\!known_failing --processes=0
	coverage html 

coverage-gcovr.xml: FORCE
	gcovr --root=. --branches --gcov-exclude='.*zlib.*|.*bzip2.*' --xml \
		--output=coverage-gcovr.xml

nosetests.xml: all
	./setup.py nosetests --with-xunit

doxygen: FORCE
	mkdir -p doc/doxygen
	sed "s/\$${VERSION}/`python ./lib/get_version.py`/" Doxyfile.in > Doxyfile
	doxygen

lib:
	cd lib && \
	$(MAKE)

test: all
	python -m nose # match the coverage command, work around bug in the setuptools
	# nose command that wipes out the build_ext config

FORCE:
