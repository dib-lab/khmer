#!/bin/bash

cov_analysis_dir=/usr/local/cov-analysis-linux64-6.6.1/bin
cov_analysis_bin=cov-build

rm -Rf .env build dist khmer/_khmermodule.so cov-int

virtualenv .env

. .env/bin/activate
pip install --quiet nose coverage pylint pep8 screed

make clean

unset coverage_pre coverage_post coverity

if [[ "${NODE_LABELS}" == *linux* ]]
then
        if type ccache >/dev/null 2>&1
        then
                echo Enabling ccache
                ccache --max-files=0 --max-size=500G
                export PATH="/usr/lib/ccache:${PATH}"
        fi
	if type gcov >/dev/null 2>&1
	then
		export CFLAGS="-pg -fprofile-arcs -ftest-coverage"
		coverage_post='--debug --inplace --libraries gcov'
	else
		echo "gcov was not found, skipping coverage check"
	fi
else
	echo "Not on a Linux node, skipping coverage check"
		
fi

if [[ "${JOB_NAME}" == khmer/* ]]
then
	if [[ -x ${cov_analysis_dir}/${cov_analysis_bin} ]]
	then
		if [[ -n "${COVERITY_TOKEN}" ]]
			#was -v COVERITY_TOKEN, but OS X bash not new enough
		then
			PATH=${PATH}:${cov_analysis_dir}
			coverity="${cov_analysis_bin} --dir cov-int"
		else
			echo "Missing coverity credentials, skipping scan"
		fi
	else
		echo "${cov_analysis_bin} does not exist in \
			${cov_analysis_dir}. Skipping coverity scan."
	fi
else
	echo "Not the main build so skipping the coverity scan"
fi

if [[ -n "${coverage_post}" ]]
then
	${coverity} python setup.py build_ext --build-temp $PWD ${coverage_post}
fi

if [[ -n "$coverity" ]]
	# was -v coverity but OS X bash not new enough
then
	tar czf khmer-cov.tgz cov-int
	curl --form project=Khmer --form token=${COVERITY_TOKEN} --form \
		email=mcrusoe@msu.edu --form file=@khmer-cov.tgz --form \
		version=`git describe --tags | sed s/v//` \
		http://scan5.coverity.com/cgi-bin/upload.py
fi
./setup.py install

./setup.py develop --build-directory .

pip install --quiet nose coverage
make coverage
make doc

	make pylint 2>&1 > pylint.out

	make pep8 2>&1 > pep8.out

	pip install --quiet -U gcovr
	make coverage-gcovr.xml

	make cppcheck-result.xml

	make doxygen 2>&1 > doxygen.out
else
	make nosetests.xml
fi
