#!/bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
for N in 0 1 2 3 4 5 6 7 8 9; do
	for M in 0 1 2 3 4 5 6 7 8 9; do
		for Z in 0 1 2 3 4 5 6 7 8 9; do
  			./ppt 0.18"$N$M" 1 1 R_0.18"$N$M"_12_R"$Z".txt 12 1.0
		done	
	done
done
octave --no-window-system makeAll.m
#./matlab -nodisplay makeAll.m
