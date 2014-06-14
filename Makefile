# make pep8 to check for basic Python code compliance
# make autopep8 to fix most pep8 errors
# make pylint to check Python code for enhanced compliance including naming
#  and documentation
# make coverage-report to check coverage of the python scripts by the tests

CPPSOURCES=$(wildcard {lib,khmer}/*.{cc,hh})
PYSOURCES=$(wildcard {khmer,scripts}/*.py)
SOURCES=$(PYSOURCES) $(CPPSOURCES) setup.py

GCOVRURL=git+https://github.com/nschum/gcovr.git@never-executed-branches
VERSION=$(shell git describe --tags --dirty | sed s/v//)

all: khmer/_khmermodule.so

khmer/_khmermodule.so: $(SOURCES)
	./setup.py build_ext --inplace

coverage-debug: $(SOURCES)
	export CFLAGS="-pg -fprofile-arcs -ftest-coverage -O0"; ./setup.py \
		build_ext --debug --inplace --libraries gcov
	touch coverage-debug

install: FORCE
	./setup.py build install

dist: dist/khmer-$(VERSION).tar.gz

dist/khmer-$(VERSION).tar.gz: $(SOURCES)
	./setup.py sdist

clean: FORCE
	cd lib && ${MAKE} clean || true
	cd tests && rm -rf khmertest_* || true
	rm -f khmer/_khmermodule.so || true
	rm khmer/*.pyc lib/*.pyc || true
	./setup.py clean --all || true
	rm coverage-debug || true

debug: FORCE
	export CFLAGS="-pg -fprofile-arcs"; python setup.py build_ext --debug \
		--inplace

doc: build/sphinx/html/index.html

build/sphinx/html/index.html: $(SOURCES) $(wildcard doc/*.txt) doc/conf.py
	pip2 install --user sphinx sphinxcontrib-autoprogram || pip2 install \
		sphinx sphinxcontrib-autoprogram
	./setup.py build_sphinx --fresh-env
	@echo ''
	@echo '--> docs in build/sphinx/html <--'
	@echo ''

pdf: build/sphinx/latex/khmer.pdf

build/sphinx/latex/khmer.pdf: $(SOURCES) doc/conf.py $(wildcard doc/*.txt)
	pip2 install --user sphinx sphinxcontrib-autoprogram || pip2 install \
		sphinx sphinxcontrib-autoprogram
	./setup.py build_sphinx --fresh-env --builder latex
	cd build/sphinx/latex && ${MAKE} all-pdf
	@echo ''
	@echo '--> pdf in build/sphinx/latex/khmer.pdf'

cppcheck-result.xml: $(CPPSOURCES)
	ls lib/*.cc khmer/_khmermodule.cc | grep -v test | cppcheck -DNDEBUG \
		-DVERSION=0.0.cppcheck -UNO_UNIQUE_RC --enable=all \
		--file-list=- -j8 --platform=unix64 --std=posix --xml \
		--xml-version=2 2> cppcheck-result.xml

cppcheck: $(CPPSOURCES)
	ls lib/*.cc khmer/_khmermodule.cc | grep -v test | cppcheck -DNDEBUG \
		-DVERSION=0.0.cppcheck -UNO_UNIQUE_RC --enable=all \
		--file-list=- -j8 --platform=unix64 --std=posix --quiet

pep8: $(PYSOURCES) $(wildcard tests/*.py)
	pip2 install --user --quiet pep8==1.5 || pip2 install --quiet pep8==1.5
	pep8 --exclude=_version.py setup.py khmer/ scripts/ tests/ || true

pep8_report.txt: $(PYSOURCES) $(wildcard tests/*.py)
	pip2 install --user --quiet pep8==1.5 || pip2 install --quiet pep8==1.5
	pep8 --exclude=_version.py setup.py khmer/ scripts/ tests/ \
		> pep8_report.txt || true

diff_pep8_report: pep8_report.txt
	pip2 install --user diff_cover || pip2 install diff_cover
	diff-quality --violations=pep8 pep8_report.txt

autopep8: $(PYSOURCES) $(wildcard tests/*.py)
	pip2 install --user autopep8 || pip2 install autopep8
	autopep8 --recursive --in-place --exclude _version.py --ignore E309 \
		setup.py khmer/ scripts/ tests/

pylint: $(PYSOURCES) $(wildcard tests/*.py)
	pip2 install --user pylint || pip2 install pylint
	pylint --msg-template="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" \
		setup.py khmer/[!_]*.py khmer/__init__.py scripts/*.py tests \
		|| true

pylint_report.txt: ${PYSOURCES} $(wildcard tests/*.py)
	pip2 install --user pylint || pip2 install pylint
	pylint --msg-template="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" \
		setup.py khmer/[!_]*.py khmer/__init__.py scripts/*.py tests \
		> pylint_report.txt || true

diff_pylint_report: pylint_report.txt
	pip2 install --user diff_cover || pip2 install diff_cover
	diff-quality --violations=pylint pylint_report.txt

# We need to get coverage to look at our scripts. Since they aren't in a
# python module we can't tell nosetests to look for them (via an import
# statement). So we run nose inside of coverage.
.coverage: $(PYSOURCES) $(wildcard tests/*.py)
	pip2 install --user --upgrade coverage || pip2 install --upgrade \
		coverage
	coverage run --branch --source=scripts,khmer --omit=khmer/_version.py \
		-m nose --with-xunit --attr=\!known_failing --processes=0

coverage.xml: .coverage
	pip2 install --user --upgrade coverage || pip2 install --upgrade \
		coverage
	coverage xml

coverage.html: htmlcov/index.html

htmlcov/index.html: .coverage
	pip2 install --user --upgrade coverage || pip2 install --upgrade \
		coverage
	coverage html
	@echo Test coverage of the Python code is now in htmlcov/index.html

coverage-report: .coverage
	pip2 install --user --upgrade coverage || pip2 install --upgrade \
		coverage
	coverage report

coverage-gcovr.xml: coverage-debug .coverage
	pip2 install --user --upgrade ${GCOVRURL}'#gcovr' || pip2 install \
		--upgrade ${GCOVRURL}'#gcovr'
	gcovr --root=. --branches --gcov-exclude='.*zlib.*|.*bzip2.*' --xml \
		--output=coverage-gcovr.xml

diff-cover: coverage-gcovr.xml coverage.xml
	pip2 install --user diff_cover || pip2 install diff_cover
	diff-cover coverage-gcovr.xml coverage.xml

diff-cover.html: coverage-gcovr.xml coverage.xml
	pip2 install --user diff_cover || pip2 install diff_cover
	diff-cover coverage-gcovr.xml coverage.xml \
		--html-report diff-cover.html

nosetests.xml: FORCE
	pip2 install --user nose || pip2 install nose
	./setup.py nosetests --with-xunit

doxygen: doc/doxygen/html/index.html

doc/doxygen/html/index.html: ${CPPSOURCES} ${PYSOURCES}
	mkdir -p doc/doxygen
	sed "s/\$${VERSION}/`python ./lib/get_version.py`/" Doxyfile.in > \
		Doxyfile
	doxygen

lib:
	cd lib && \
	$(MAKE)

test:
	pip2 install --user nose || pip2 install nose
	./setup.py nosetests

sloccount.sc: ${CPPSOURCES} ${PYSOURCES} $(wildcard tests/*.py) Makefile
	sloccount --duplicates --wide --details lib khmer scripts tests \
		setup.py Makefile > sloccount.sc

sloccount: 
	sloccount lib khmer scripts tests setup.py Makefile

coverity-build:
	if [[ -x ${cov_analysis_dir}/bin/cov-build ]]; \
	then \
		export PATH=${PATH}:${cov_analysis_dir}/bin; \
		cov-build --dir cov-int --c-coverage gcov --disable-gcov-arg-injection make coverage-debug; \
		cov-capture --dir cov-int --c-coverage gcov python -m nose --attr '!known_failing' ; \
		cov-import-scm --dir cov-int --scm git 2>/dev/null; \
	else echo 'bin/cov-build does not exist in $$cov_analysis_dir: '\
		'${cov_analysis_dir}. Skipping coverity scan.'; \
	fi

coverity-upload: cov-int
	if [[ -n "${COVERITY_TOKEN}" ]]; \
	then \
		tar czf khmer-cov.tgz cov-int; \
		curl --form project=ged-lab/khmer \
			--form token=${COVERITY_TOKEN} --form \
			email=mcrusoe@msu.edu --form file=@khmer-cov.tgz \
			--form version=${VERSION} \
			http://scan5.coverity.com/cgi-bin/upload.py; \
	else echo 'Missing coverity credentials in $$COVERITY_TOKEN,'\
		'skipping scan'; \
	fi

coverity-clean-configuration:
	rm -f ${cov_analysis_dir}/config/coverity_config.xml

coverity-configure:
	if [[ -x ${cov_analysis_dir}/bin/cov-configure ]]; \
	then \
		export PATH=${PATH}:${cov_analysis_dir}/bin; \
		for compiler in /usr/bin/gcc-4.8 /usr/bin/x86_64-linux-gnu-gcc; do \
       			cov-configure --comptype gcc --compiler $${compiler}; \
		done; \
	else echo 'bin/cov-configure does not exist in $$cov_analysis_dir: '\
		'${cov_analysis_dir}. Skipping coverity configuration.'; \
	fi

FORCE:
