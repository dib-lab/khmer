#!/bin/bash
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
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

make clean

rm -Rf .env dist cov-int

if [ -z "${PYTHON_EXECUTABLE}" ]; then
    if type python2> /dev/null 2>&1
    then
        PYTHON_EXECUTABLE=$(which python2)
    else
        PYTHON_EXECUTABLE=$(which python)
    fi
fi
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
	make coverage-gcovr.xml coverage.xml TESTATTR='"not known_failing and not huge"'
	./setup.py install
else
	echo "gcov was not found (or we are on OSX), skipping coverage check"
	./setup.py install
	./setup.py develop
	make pytests.xml
fi

if type cppcheck >/dev/null 2>&1
then
	make cppcheck-result.xml SHELL=bash
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

unset CFLAGS
unset LDFLAGS
unset CPPFLAGS
unset CXXFLAGS

# Don't do lib too, as we already compile as part of libtest
make libtest

# Upload code coverage to codecov.io
pip install codecov
codecov -X pycov search gcov -f coverage.xml coverage-gcovr.xml &> /dev/null
