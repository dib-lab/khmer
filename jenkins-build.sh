#!/bin/bash

make clean

rm -Rf .env dist cov-int

virtualenv .env

. .env/bin/activate
pip install --quiet nose coverage pylint pep8==1.5 screed

if type ccache >/dev/null 2>&1
then
        echo Enabling ccache
        ccache --max-files=0 --max-size=500G
        export PATH="/usr/lib/ccache:${PATH}"
fi
if [[ "${NODE_LABELS}" == *osx* ]]
then
	export ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future
fi

if type gcov >/dev/null 2>&1 && [[ "${NODE_LABELS}" != *osx* ]]
then
	export CFLAGS="-pg -fprofile-arcs -ftest-coverage"
	python setup.py build_ext --build-temp $PWD --debug --inplace \
		--libraries gcov
	make coverage-gcovr.xml coverage.xml
	./setup.py install
else
	echo "gcov was not found (or we are on OSX), skipping coverage check"
	./setup.py install
	make nosetests.xml
fi

if type cppcheck >/dev/null 2>&1
then
	make cppcheck-result.xml
fi
if type doxygen >/dev/null 2>&1
then
	make doxygen 2>&1 > doxygen.out
fi

if type hg >/dev/null 2>&1
then
	rm -Rf sphinx-contrib
	#hg clone http://bitbucket.org/mcrusoe/sphinx-contrib
	hg clone http://athyra.ged.msu.edu/~mcrusoe/sphinx-contrib
	pip install --upgrade sphinx-contrib/autoprogram/
	make doc
fi
make pylint 2>&1 > pylint.out
make pep8 2>&1 > pep8.out

if type sloccount >/dev/null 2>&1
then
	make sloccount.sc
fi
