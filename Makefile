# make pep8 to check for basic Python code compliance
# make autopep8 to fix most pep8 errors
# make pylint to check Python code for enhanced compliance including naming
#  and documentation
# make coverage-report to check coverage of the python scripts by the tests

CPPSOURCES=$(wildcard lib/*.cc lib/*.hh khmer/_khmermodule.cc)
PYSOURCES=$(wildcard khmer/*.py scripts/*.py)
SOURCES=$(PYSOURCES) $(CPPSOURCES) setup.py
DEVPKGS=sphinxcontrib-autoprogram pep8==1.5.7 diff_cover \
autopep8 pylint coverage gcovr nose pep257 future screed

GCOVRURL=git+https://github.com/nschum/gcovr.git@never-executed-branches
VERSION=$(shell git describe --tags --dirty | sed s/v//)
CPPCHECK=ls lib/*.cc khmer/_khmermodule.cc | grep -v test | cppcheck -DNDEBUG \
	 -DVERSION=0.0.cppcheck -UNO_UNIQUE_RC --enable=all \
	 --file-list=- --platform=unix64 --std=c++03 --inline-suppr \
	 --quiet -Ilib -Ithird-party/bzip2 -Ithird-party/zlib \
	 -Ithird-party/smhasher

UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
	TESTATTR='!known_failing,!jenkins'
else
	TESTATTR='!known_failing,!jenkins,!linux'
endif

## all         : default task; compile C++ code, build shared object library
all: sharedobj

## help        : print this help message and exit
help: Makefile
	@sed -n 's/^##//p' $<

## install-dep : install most of the development dependencies via pip
install-dep: install-dependencies

install-dependencies:
	pip install --upgrade $(DEVPKGS)

## sharedobj   : build khmer shared object file
sharedobj: khmer/_khmermodule.so

khmer/_khmermodule.so: $(CPPSOURCES)
	./setup.py build_ext --inplace

coverage-debug: $(CPPSOURCES)
	export CFLAGS="-pg -fprofile-arcs -ftest-coverage -O0"; ./setup.py \
		build_ext --debug --inplace --libraries gcov
	touch coverage-debug

## install     : install the khmer module and scripts
install: FORCE
	./setup.py build install

## dist        : create a module package for distribution
dist: dist/khmer-$(VERSION).tar.gz

dist/khmer-$(VERSION).tar.gz: $(SOURCES)
	./setup.py sdist

## clean       : clean up all temporary / machine-generated files
clean: FORCE
	cd lib && ${MAKE} clean || true
	cd tests && rm -rf khmertest_* || true
	rm -f khmer/_khmermodule.so
	rm -f khmer/*.pyc lib/*.pyc
	./setup.py clean --all || true
	rm -f coverage-debug
	rm -Rf .coverage
	rm -f diff-cover.html

debug: FORCE
	export CFLAGS="-pg -fprofile-arcs"; python setup.py build_ext --debug \
		--inplace

## doc         : render documentation in HTML
doc: build/sphinx/html/index.html

build/sphinx/html/index.html: $(SOURCES) $(wildcard doc/*.txt) doc/conf.py all
	./setup.py build_sphinx --fresh-env
	@echo ''
	@echo '--> docs in build/sphinx/html <--'
	@echo ''

## pdf         : render documentation as a PDF file
pdf: build/sphinx/latex/khmer.pdf

build/sphinx/latex/khmer.pdf: $(SOURCES) doc/conf.py $(wildcard doc/*.txt)
	./setup.py build_sphinx --fresh-env --builder latex
	cd build/sphinx/latex && ${MAKE} all-pdf
	@echo ''
	@echo '--> pdf in build/sphinx/latex/khmer.pdf'

cppcheck-result.xml: $(CPPSOURCES)
	${CPPCHECK} --xml-version=2 2> cppcheck-result.xml

## cppcheck    : run static analysis on C++ code
cppcheck: $(CPPSOURCES)
	${CPPCHECK}

## pep8        : check Python code style
pep8: $(PYSOURCES) $(wildcard tests/*.py)
	pep8 --exclude=_version.py  --show-source --show-pep8 setup.py khmer/ \
		scripts/ tests/ oxli/ || true

pep8_report.txt: $(PYSOURCES) $(wildcard tests/*.py)
	pep8 --exclude=_version.py setup.py khmer/ scripts/ tests/ oxli/ \
		> pep8_report.txt || true

diff_pep8_report: pep8_report.txt
	diff-quality --violations=pep8 pep8_report.txt

## pep257      : check Python code style
pep257: $(PYSOURCES) $(wildcard tests/*.py)
	pep257 --ignore=D100,D101,D102,D103 \
		setup.py khmer/ scripts/ tests/ || true

pep257_report.txt: $(PYSOURCES) $(wildcard tests/*.py)
	pep257 setup.py khmer/ scripts/ tests/ \
		> pep257_report.txt 2>&1 || true

diff_pep257_report: pep257_report.txt
	diff-quality --violations=pep8 pep257_report.txt

## astyle      : fix most C++ code indentation and formatting
astyle: $(CPPSOURCES)
	astyle -A10 --max-code-length=80 $(CPPSOURCES)

## autopep8    : fix most Python code indentation and formatting
autopep8: $(PYSOURCES) $(wildcard tests/*.py)
	autopep8 --recursive --in-place --exclude _version.py --ignore E309 \
		setup.py khmer/*.py scripts/*.py tests/*.py oxli/*.py

# A command to automatically run astyle and autopep8 on appropriate files
## format      : check/fix all code indentation and formatting (runs astyle and autopep8)
format: astyle autopep8
	# Do nothing

## pylint      : run static code analysis on Python code
pylint: $(PYSOURCES) $(wildcard tests/*.py)
	pylint --msg-template="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" \
		setup.py khmer/[!_]*.py khmer/__init__.py scripts/*.py tests \
		oxli/*.py || true

pylint_report.txt: ${PYSOURCES} $(wildcard tests/*.py)
	pylint --msg-template="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" \
		setup.py khmer/[!_]*.py khmer/__init__.py scripts/*.py tests \
		sandbox/*.py oxli/*.py > pylint_report.txt || true

diff_pylint_report: pylint_report.txt
	diff-quality --violations=pylint pylint_report.txt

# We need to get coverage to look at our scripts. Since they aren't in a
# python module we can't tell nosetests to look for them (via an import
# statement). So we run nose inside of coverage.
.coverage: $(PYSOURCES) $(wildcard tests/*.py) khmer/_khmermodule.so
	coverage run --branch --source=scripts,khmer,oxli --omit=khmer/_version.py \
		-m nose --with-xunit --attr=\!known_failing --processes=0

coverage.xml: .coverage
	coverage xml

coverage.html: htmlcov/index.html

htmlcov/index.html: .coverage
	coverage html
	@echo Test coverage of the Python code is now in htmlcov/index.html

coverage-report: .coverage
	coverage report

coverage-gcovr.xml: coverage-debug .coverage
	gcovr --root=. --branches --output=coverage-gcovr.xml --xml \
          --gcov-exclude='.*zlib.*|.*bzip2.*|.*smhasher.*|.*seqan.*'

diff-cover: coverage-gcovr.xml coverage.xml
	diff-cover coverage-gcovr.xml coverage.xml

diff-cover.html: coverage-gcovr.xml coverage.xml
	diff-cover coverage-gcovr.xml coverage.xml \
		--html-report diff-cover.html

nosetests.xml: FORCE
	./setup.py nosetests --with-xunit --attr ${TESTATTR}

## doxygen     : generate documentation of the C++ and Python code
doxygen: doc/doxygen/html/index.html

doc/doxygen/html/index.html: ${CPPSOURCES} ${PYSOURCES}
	mkdir -p doc/doxygen
	sed "s/\$${VERSION}/`python ./lib/get_version.py`/" Doxyfile.in > \
		Doxyfile
	doxygen

lib: FORCE
	cd lib && \
	$(MAKE)

# Runs a test of ./lib
libtest: FORCE
	rm -rf install_target
	mkdir -p install_target
	cd lib && \
	 $(MAKE) clean && \
	 $(MAKE) all && \
	 $(MAKE) install PREFIX=../install_target
	test -d install_target/include
	test -f install_target/include/khmer.hh
	test -d install_target/lib
	test -f install_target/lib/libkhmer.a
	$(CXX) -o install_target/test-prog-static -I install_target/include \
		lib/test-compile.cc install_target/lib/libkhmer.a
	$(CXX) -o install_target/test-prog-dynamic -I install_target/include \
		-L install_target/lib lib/test-compile.cc -lkhmer
	rm -rf install_target

## test        : run the khmer test suite
test: FORCE
	./setup.py develop
	./setup.py nosetests --attr ${TESTATTR}

sloccount.sc: ${CPPSOURCES} ${PYSOURCES} $(wildcard tests/*.py) Makefile
	sloccount --duplicates --wide --details lib khmer scripts tests \
		setup.py Makefile > sloccount.sc

## sloccount   : count lines of code
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

compile_commands.json: clean
	export PATH=$(shell echo $$PATH | sed 's=/usr/lib/ccache:==g') ; \
		bear -- ./setup.py build_ext

convert-release-notes:
	for file in doc/release-notes/*.md; do \
		pandoc --from=markdown --to=rst $${file} > $${file%%.md}.rst; \
		done

FORCE:
