#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
cd /root/screed
git pull origin master
python setup.py install

cd /root/khmer
make clean
git pull origin refactor
make clean test
