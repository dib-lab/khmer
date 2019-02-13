# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org

# `make help` for a summary of useful targets

# `SHELL=bash` Will break Titus's laptop, so don't use BASH-isms like
# `[[` conditional expressions.
#
PREFIX=/usr/local
CPPSOURCES=$(wildcard src/oxli/*.cc include/oxli/*.hh src/khmer/_cpy_*.cc include/khmer/_cpy_*.hh) setup.py
CYSOURCES=$(wildcard khmer/_oxli/*.pxd khmer/_oxli/*.pyx)
PYSOURCES=$(filter-out khmer/_version.py, \
	  $(wildcard khmer/*.py scripts/*.py oxli/*.py) )
SOURCES=$(PYSOURCES) $(CPPSOURCES) $(CYSOURCES) setup.py

DEVPKGS=pep8==1.6.2 diff_cover autopep8 pylint coverage gcovr pytest pydocstyle

GCOVRURL=git+https://github.com/nschum/gcovr.git@never-executed-branches

VERSION=$(shell ./setup.py version | grep Version | awk '{print $$4}' \
	| sed 's/+/-/')

# The following four variables are only used by cppcheck. If you want to
# change how things are compiled edit `setup.cfg` or `setup.py`.
DEFINES += -DNDEBUG -DVERSION=$(VERSION) -DSEQAN_HAS_BZIP2=1 \
	   -DSEQAN_HAS_ZLIB=1 -UNO_UNIQUE_RC

INCLUDESTRING=$(shell gcc -E -x c++ - -v < /dev/null 2>&1 >/dev/null \
	    | grep '^ /' | grep -v cc1plus)
INCLUDEOPTS=$(shell gcc -E -x c++ - -v < /dev/null 2>&1 >/dev/null \
	    | grep '^ /' | grep -v cc1plus | awk '{print "-I" $$1 " "}')

PYINCLUDE=$(shell python -c "import sysconfig; \
            flags = ['-I' + sysconfig.get_path('include'), \
            '-I' + sysconfig.get_path('platinclude')]; print(' '.join(flags))")

CPPCHECK_SOURCES=$(filter-out lib/test%, $(wildcard lib/*.cc khmer/_khmer.cc) )
CPPCHECK=cppcheck --enable=all \
	 --error-exitcode=1 \
	 --suppress='*:/Library/*' \
	 --suppress='*:*/include/python*/Python.h' \
	 --suppress='*:/usr/*' --platform=unix64 \
	 --std=c++11 --inline-suppr -Ilib -Ithird-party/bzip2 \
	 -Ithird-party/zlib -Ithird-party/smhasher -Ithird-party/rollinghash \
	 $(DEFINES) $(INCLUDEOPTS) $(PYINCLUDE) $(CPPCHECK_SOURCES) --quiet

UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
	TESTATTR ?= 'not known_failing and not jenkins and not huge'
else
	TESTATTR ?= 'not known_failing and not jenkins and not huge and not linux'
endif

MODEXT=$(shell python -c \
       "import sysconfig;print(sysconfig.get_config_var('SO'))")
EXTENSION_MODULE = khmer/_khmer$(MODEXT)
CY_MODULES = $($(wildcard khmer/_oxli/*.pyx): .pyx=.$(MODEXT))

PYLINT_TEMPLATE="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}"

## all         : default task; compile C++ code, build shared object library
all: sharedobj

## help        : print this help message and exit
help: Makefile
	@sed -n 's/^##//p' $<

## install-dep : install most of the development dependencies via pip
install-dep: install-dependencies

install-dependencies:
	pip install $(DEVPKGS)
	pip install --requirement doc/requirements.txt

## sharedobj   : build khmer shared object file
sharedobj: $(EXTENSION_MODULE)

$(EXTENSION_MODULE): $(CPPSOURCES) $(CYSOURCES)
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
	cd src/oxli && $(MAKE) clean || true
	cd tests && rm -rf khmertest_* || true
	rm -f pytests.xml
	cd third-party/cqf && make clean || true
	rm -f $(EXTENSION_MODULE)
	rm -f khmer/*.pyc scripts/*.pyc tests/*.pyc oxli/*.pyc \
		sandbox/*.pyc khmer/__pycache__/* sandbox/__pycache__/* \
		khmer/_oxli/*.so
	./setup.py clean --all || true
	rm -f coverage-debug
	rm -Rf .coverage coverage-gcovr.xml coverage.xml
	rm -f diff-cover.html
	rm -Rf build dist
	rm -rf __pycache__/ khmer.egg-info/
	@find ./ -type d -name __pycache__ -exec rm -rf {} +
	@find ./khmer/ -type f -name *$(MODEXT) -exec rm -f {} +
	-rm -f *.gcov

debug: FORCE
	export CFLAGS="-pg -fprofile-arcs -D_GLIBCXX_DEBUG_PEDANTIC \
		-D_GLIBCXX_DEBUG -DDEBUG_ASSEMBLY=1 -DDEBUG_FILTERS=1"; python setup.py build_ext --debug \
		--inplace

## doc         : render documentation in HTML
doc: build/sphinx/html/index.html

build/sphinx/html/index.html: $(SOURCES) $(wildcard doc/*.rst) doc/conf.py all
	./setup.py build_sphinx --fresh-env
	@echo ''
	@echo '--> docs in build/sphinx/html <--'
	@echo ''

## pdf         : render documentation as a PDF file
# packages needed include: texlive-latex-base texlive-latex-recommended
# texlive-fonts-recommended texlive-latex-extra
pdf: build/sphinx/latex/khmer.pdf

build/sphinx/latex/khmer.pdf: $(SOURCES) doc/conf.py $(wildcard doc/*.rst) \
	$(wildcard doc/user/*.rst) $(wildcard doc/dev/*.rst) sharedobj
	./setup.py build_sphinx --fresh-env --builder latex
	cd build/sphinx/latex && $(MAKE) all-pdf
	@echo ''
	@echo '--> pdf in build/sphinx/latex/khmer.pdf'

cppcheck-result.xml: $(CPPSOURCES)
	$(CPPCHECK) --xml-version=2 2> cppcheck-result.xml

## cppcheck    : run static analysis on C++ code
cppcheck: FORCE
	@$(CPPCHECK)

cppcheck-long: FORCE
	@$(CPPCHECK) -Ithird-party/seqan/core/include

## pep8        : check Python code style
pep8: $(PYSOURCES) $(wildcard tests/*.py)
	pep8 setup.py khmer/*.py scripts/*.py tests/*.py oxli/*.py examples/python-api/*.py

pep8_report.txt: $(PYSOURCES) $(wildcard tests/*.py)
	pep8 setup.py khmer/ scripts/ tests/ oxli/ \
		> pep8_report.txt || true

diff_pep8_report: pep8_report.txt
	diff-quality --violations=pep8 pep8_report.txt

## pydocstyle  : check Python doc strings
pydocstyle: $(PYSOURCES) $(wildcard tests/*.py)
	pydocstyle --ignore=D100,D101,D102,D103,D203 --match='(?!_version).*\.py' \
		setup.py khmer/ scripts/ oxli/ || true

pydocstyle_report.txt: $(PYSOURCES) $(wildcard tests/*.py)
	pydocstyle setup.py khmer/ scripts/ tests/ oxli/ \
		> pydocstyle_report.txt 2>&1 || true

diff_pydocstyle_report: pydocstyle_report.txt
	diff-quality --violations=pep8 pydocstyle_report.txt

## astyle      : fix most C++ code indentation and formatting
astyle: $(CPPSOURCES)
	astyle -A10 --max-code-length=80 $(filter-out setup.py,$(CPPSOURCES))

## autopep8    : fix most Python code indentation and formatting
autopep8: $(PYSOURCES) $(wildcard tests/*.py)
	autopep8 --recursive --in-place --exclude _version.py --ignore E309 \
		setup.py khmer/*.py scripts/*.py tests/*.py oxli/*.py

## format      : check/fix all code formatting (astyle and autopep8)
format: astyle autopep8
	# Do nothing

## pylint      : run static code analysis on Python code
pylint: $(PYSOURCES) $(wildcard tests/*.py)
	pylint --msg-template=$(PYLINT_TEMPLATE) \
		setup.py $(PYSOURCES) tests/*.py || true

pylint_report.txt: $(PYSOURCES) $(wildcard tests/*.py) $(wildcard sandbox/*.py)
	pylint --msg-template=$(PYLINT_TEMPLATE) \
		$(PYSOURCES) tests sandbox > pylint_report.txt || true

diff_pylint_report: pylint_report.txt
	diff-quality --violations=pylint pylint_report.txt

# We need to get coverage to look at our scripts. Since they aren't in a
# python module we can't tell pytest to look for them (via an import
# statement). So we run pytest inside of coverage.
.coverage: $(PYSOURCES) $(wildcard tests/*.py) $(EXTENSION_MODULE)
	./setup.py develop
	coverage run --branch --source=scripts,khmer,oxli \
		--omit=khmer/_version.py -m pytest --junitxml=pytests.xml \
		-m $(TESTATTR)

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
          --gcov-exclude='.*zlib.*|.*bzip2.*|.*smhasher.*|.*seqan.*' \
	  --exclude-unreachable-branches

diff-cover: coverage-gcovr.xml coverage.xml
	diff-cover coverage-gcovr.xml coverage.xml

diff-cover.html: coverage-gcovr.xml coverage.xml
	diff-cover coverage-gcovr.xml coverage.xml \
		--html-report diff-cover.html

pytests.xml: FORCE
	py.test --junitxml=$@ -m ${TESTATTR}

## doxygen     : generate documentation of the C++ and Python code
# helpful packages: doxygen graphviz
# ignore warning re: _formulas.aux
doxygen: doc/doxygen/html/index.html

doc/doxygen/html/index.html: $(CPPSOURCES) $(PYSOURCES)
	mkdir -p doc/doxygen
	sed "s=\$${VERSION}=$(VERSION)=" Doxyfile.in | \
		sed "s=\$${INCLUDES}=$(INCLUDESTRING)=" > Doxyfile
	doxygen

liboxli: FORCE
	cd src/oxli && \
	$(MAKE)

install-liboxli: liboxli
	cd src/oxli && $(MAKE) install PREFIX=$(PREFIX)
	mkdir -p $(PREFIX)/include/khmer
	cp -r include/khmer/_cpy_*.hh $(PREFIX)/include/khmer/

# Runs a test of liboxli
libtest: FORCE
	rm -rf install_target
	mkdir -p install_target
	cd src/oxli && \
	 $(MAKE) clean && \
	 $(MAKE) all && \
	 $(MAKE) install PREFIX=../install_target
	test -d install_target/include
	test -f install_target/include/oxli/oxli.hh
	test -d install_target/lib
	test -f install_target/lib/liboxli.a
	$(CXX) -std=c++11 -o install_target/test-prog-static \
		-I install_target/include src/oxli/test-compile.cc \
		install_target/lib/liboxli.a
	$(CXX) -std=c++11 -o install_target/test-prog-dynamic \
		-I install_target/include -L install_target/lib \
		src/oxli/test-compile.cc -loxli
	rm -rf install_target

## test        : run the khmer test suite
test: FORCE
	./setup.py develop
	py.test -m ${TESTATTR}

sloccount.sc: $(CPPSOURCES) $(PYSOURCES) $(wildcard tests/*.py) Makefile
	sloccount --duplicates --wide --details include src khmer scripts tests \
		setup.py Makefile > sloccount.sc

## sloccount   : count lines of code
sloccount:
	sloccount src include khmer scripts tests setup.py Makefile

coverity-build:
	if [ -x "${cov_analysis_dir}/bin/cov-build" ]; \
	then \
		export PATH=${PATH}:${cov_analysis_dir}/bin; \
		cov-build --dir cov-int --c-coverage gcov \
			--disable-gcov-arg-injection make coverage-debug; \
		cov-capture --dir cov-int --c-coverage gcov python -m pytest \
			-m $(TESTATTR) ; \
		cov-import-scm --dir cov-int --scm git 2>/dev/null; \
	else echo 'bin/cov-build does not exist in $$cov_analysis_dir: '\
		'${cov_analysis_dir}. Skipping coverity scan.'; \
	fi

coverity-upload: cov-int
	if [ -n "${COVERITY_TOKEN}" ]; \
	then \
		tar czf khmer-cov.tgz cov-int; \
		curl --form token=${COVERITY_TOKEN} --form \
			email=mcrusoe@msu.edu --form file=@khmer-cov.tgz \
			--form version=$(VERSION) \
			https://scan.coverity.com/builds?project=ged-lab%2Fkhmer ; \
	else echo 'Missing coverity credentials in $$COVERITY_TOKEN,'\
		'skipping scan'; \
	fi

coverity-clean-configuration:
	rm -f ${cov_analysis_dir}/config/coverity_config.xml

coverity-configure:
	if [[ -x ${cov_analysis_dir}/bin/cov-configure ]]; \
	then \
		export PATH=${PATH}:${cov_analysis_dir}/bin; \
		for compiler in \
			/usr/bin/gcc-4.8 /usr/bin/x86_64-linux-gnu-gcc; do \
			cov-configure --comptype gcc \
					--compiler $${compiler}; \
		done; \
	else echo 'bin/cov-configure does not exist in $$cov_analysis_dir: '\
		'${cov_analysis_dir}. Skipping coverity configuration.'; \
	fi

# may need to `sudo apt-get install bear`
compile_commands.json: clean
	export PATH=$(shell echo $$PATH | sed 's=/usr/lib/ccache:==g') ; \
		bear ./setup.py build_ext

convert-release-notes:
	for file in doc/release-notes/*.md; do \
		pandoc --from=markdown --to=rst $${file} > $${file%%.md}.rst; \
		done

list-author-emails:
	@echo 'name,E-Mail Address'
	@echo 'Daniel Standage,daniel.standage@gmail.com'
	@git log --format='%aN,%aE' | sort -u | grep -v -F -f author-skips.txt
	@echo 'C. Titus Brown,ctbrown@ucdavis.edu'

list-citation:
	git log --format='%aN,%aE' | sort -u | grep -v -F -f author-skips.txt > authors.csv
	python sort-authors-list.py

## cpp-demos   : run programs demonstrating access to the (unstable) C++ API
cpp-demos: sharedobj
	cd examples/c++-api/ && make all run

## py-demos    : run programs demonstrating access to the Python API
py-demos: sharedobj
	python examples/python-api/exact-counting.py
	python examples/python-api/bloom.py
	python examples/python-api/consume.py examples/c++-api/reads.fastq
	python examples/python-api/mask.py

COMMIT ?= $(shell git rev-parse HEAD)
SLUG ?= $(TRAVIS_PULL_REQUEST_SLUG)
docker-container:
	cd docker && docker build --build-arg branch=$(COMMIT) --build-arg slug=$(SLUG) .

FORCE:

# Use this to print the value of a Makefile variable
# Example `make print-VERSION`
# From https://www.cmcrossroads.com/article/printing-value-makefile-variable
print-%  : ; @echo $* = $($*)
