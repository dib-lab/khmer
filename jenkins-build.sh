#!/bin/bash

make clean

rm -Rf .env dist cov-int

#if type python2> /dev/null 2>&1
#then
#    PYTHON_EXECUTABLE=$(which python2)
#else
#    PYTHON_EXECUTABLE=$(which python)
#fi
virtualenv -p ${PYTHON_EXECUTABLE} .env

. .env/bin/activate
pip install setuptools==3.4.1
make install-dependencies

if type ccache >/dev/null 2>&1
then
        echo Enabling ccache
        ccache --max-files=0 --max-size=500G
        export PATH="/usr/lib/ccache:${PATH}"
fi
if [[ "${NODE_LABELS}" == *osx* ]]
then
	export ARCHFLAGS=-Wno-error=unused-command-line-argument
fi

if type gcov >/dev/null 2>&1 && [[ "${NODE_LABELS}" != *osx* ]]
then
	export CFLAGS="-pg -fprofile-arcs -ftest-coverage"
	python setup.py build_ext --build-temp $PWD --debug --inplace \
		--libraries gcov develop
	make coverage-gcovr.xml coverage.xml
	./setup.py install
else
	echo "gcov was not found (or we are on OSX), skipping coverage check"
	./setup.py install
	./setup.py develop
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
	#hg clone http://athyra.ged.msu.edu/~mcrusoe/sphinx-contrib
	#pip install --upgrade sphinx-contrib/autoprogram/
	#pip install -r doc/requirements.txt # now covered by make install-dep
	make doc
fi
make pylint 2>&1 > pylint.out
make pep8 2>&1 > pep8.out

if type sloccount >/dev/null 2>&1
then
	make sloccount.sc
fi

# takes too long to run on every build
#bash -ex -c 'cd examples/stamps/; ./do.sh' || { echo examples/stamps/do.sh no longer runs; /bin/false; }

make lib
make libtest
