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
	cd lib && ${MAKE} clean
	cd tests && rm -rf khmertest_*
	rm -f khmer/_khmermodule.so
	./setup.py clean --all

debug:
	export CFLAGS="-pg -fprofile-arcs"; python setup.py build_ext --debug \
		--inplace

doc: FORCE
	pip install --user sphinx || pip install sphinx
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
	pip install --user --quiet pep8==1.5 || pip install --quiet pep8==1.5
	pep8 --exclude=_version.py setup.py khmer/ scripts/ tests/ || true

autopep8: FORCE
	pip install --user autopep8 || pip install autopep8
	autopep8 --recursive --in-place --exclude _version.py --ignore E309 setup.py \
		khmer/ scripts/ tests/

pylint: FORCE
	pip install --user pylint || pip install pylint
	pylint -f parseable setup.py khmer/[!_]*.py khmer/__init__.py scripts/*.py \
		tests || true

# We need to get coverage to look at our scripts. Since they aren't in a
# python module we can't tell nosetests to look for them (via an import
# statement). So we run nose inside of coverage.
.coverage: FORCE
	pip install --user coverage || pip install coverage
	coverage run --branch --source=scripts,khmer -m nose --with-xunit \
		--attr=\!known_failing --processes=0

coverage.xml: .coverage
	coverage xml

coverage.html: .coverage
	coverage html

coverage-gcovr.xml: FORCE
	pip install --user gcovr || pip install gcovr
	gcovr --root=. --branches --gcov-exclude='.*zlib.*|.*bzip2.*' --xml \
		--output=coverage-gcovr.xml

nosetests.xml: all
	pip install --user nose || pip install nose
	./setup.py nosetests --with-xunit

doxygen: FORCE
	mkdir -p doc/doxygen
	sed "s/\$${VERSION}/`python ./lib/get_version.py`/" Doxyfile.in > Doxyfile
	doxygen

lib:
	cd lib && \
	$(MAKE)

test: all
	pip install --user nose || pip install nose
	./setup.py nosetests

FORCE:
