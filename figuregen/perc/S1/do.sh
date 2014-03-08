#!/bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
for N in 0 1 2 3 4 5 6 7 8 9; do
	for M in 0 1 2 3 4 5 6 7 8 9; do
		for Z in 4 5 6 7 8 9 10 11 12; do
#		for Z in 4 5 6 7; do
			# to run more repeats (use 100)
  			./ppt 0."$N$M" 1 100 R_0."$N$M"_"$Z".txt $Z 1.0
  			#./ppt 0."$N$M" 1 2 R_0."$N$M"_"$Z".txt $Z 1.0
		done	
	done
done
octave --no-window-system makeAll.m
#matlab -nodisplay makeAll.m
