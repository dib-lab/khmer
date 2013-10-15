#!/bin/bash

cov_analysis_dir=/usr/local/cov-analysis-linux64-6.5.1/bin
cov_analysis_bin=cov-build

rm -Rf .env build dist khmer/_khmermodule.so

virtualenv .env

. .env/bin/activate

python setup.py clean --all

unset coverage_pre coverage_post coverity

if [[ "${NODE_LABELS}" == *linux* ]]
then
	if type gcov >/dev/null 2>&1
	then
		coverage_pre='CFLAGS="-pg -fprofile-arcs -ftest-coverage"'
		coverage_post='--debug --inplace --libraries gcov'
	else
		echo "gcov was not found, skipping coverage check"
	fi
else
	echo "Not on a Linux node, skipping coverage check"
		
fi

if [[ "${JOB_NAME}" == "khmer-multi" ]]
then
	if [[ -x ${cov_analysis_dir}/${cov_analysis_bin} ]]
	then
		if [[ -v COVERITY_TOKEN ]]
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

${coverage_pre} ${coverity} python setup.py build_ext ${coverage_post}

if [[ -v coverity ]]
then
	tar czf khmer-cov.tgz cov-int
	curl --form project=Khmer --form token=${COVERITY_TOKEN} --form \
		email=mcrusoe@msu.edu --form file=@khmer-cov.tgz --form \
		version=0.0.0.${BUILD_TAG} \
		http://scan5.coverity.com/cgi-bin/upload.py
fi


pip install --quiet nosexcover
python setup.py nosetests --with-xcoverage --with-xunit --cover-package=khmer \
	--cover-erase --attr=\!known_failing

make doc

pip install --quiet pylint
pylint -f parseable doc/*.py figuregen/*.py novelty/*.py khmer/*.py sandbox/*.py \
       	scripts/*.py tests khmer | tee ../pylint.out

if [[ -v coverage_pre ]]
then
	pip install -U gcovr
	gcovr -r $PWD --xml > coverage-gcovr.xml

	cppcheck --std=posix --platform=unix64 -j8 --enable=all -I lib/ -i lib/zlib/ \
		-i lib/bzip2/ -DVALIDATE_PARTITIONS --xml lib 2> cppcheck-result.xml

	mkdir -p doc/doxygen
	doxygen
fi
