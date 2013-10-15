#!/bin/bash

set -x #print commands as they are run

rm -Rf .env build dist khmer/_khmermodule.so
virtualenv .env

. .env/bin/activate

python setup.py clean --all

#PATH=$PATH:/usr/local/cov-analysis-linux64-6.5.1/bin
#CFLAGS="-pg -fprofile-arcs -ftest-coverage" cov-build --dir cov-int python setup.py build_ext --debug 
#tar czf khmer-cov.tgz cov-int
#curl --form project=Khmer --form token=${COVERITY_TOKEN} --form email=mcrusoe@msu.edu --form file=@khmer-cov.tgz --form version=0.4.0.0.BUILD_TAG http://scan5.coverity.com/cgi-bin/upload.py

if [[ "${NODE_LABELS}" == *linux* ]]
then
	CFLAGS="-pg -fprofile-arcs -ftest-coverage" python setup.py build_ext \
		--debug --inplace --libraries gcov
else
	python setup.py build_ext
fi

pip install --quiet nosexcover
python setup.py nosetests --with-xcoverage --with-xunit --cover-package=khmer \
	--cover-erase --attr=\!known_failing

make doc

pip install --quiet pylint
pylint -f parseable doc/*.py figuregen/*.py novelty/*.py khmer/*.py sandbox/*.py \
       	scripts/*.py tests khmer | tee ../pylint.out

if [[ "${NODE_LABELS}" == *linux* ]]
then
	pip install -U gcovr
	gcovr -r $PWD --xml > coverage-gcovr.xml

	cppcheck --std=posix --platform=unix64 -j8 --enable=all -I lib/ -i lib/zlib/ \
		-i lib/bzip2/ -DVALIDATE_PARTITIONS --xml lib 2> cppcheck-result.xml

	mkdir -p doc/doxygen
	doxygen
fi
