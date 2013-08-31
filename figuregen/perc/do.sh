#!/bin/sh
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

cd src
make

cp ppt ../scripts
cp ppt ../S1

cd ../scripts
./do.sh

cd ../S1
./do.sh
