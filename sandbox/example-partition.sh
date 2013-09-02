#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
READ_FILE=data/25k.fa

export PYTHONPATH=python/
env/bin/python ./scripts/do-th-subset-save.py $READ_FILE
