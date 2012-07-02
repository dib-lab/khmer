#!/bin/bash
for N in 0 1 2 3 4 5 6 7 8 9; do
	for M in 0 1 2 3 4 5 6 7 8 9; do
#		for Z in 4 5 6 7 8 9 10 11 12; do
		for Z in 4 5 6 7; do
			# to run more repeats (use 100)
  			# ./PPT_graph 0."$N$M" 1 100 R_0."$N$M"_K"$Z".txt $Z 1.0
  			./PPT_graph 0."$N$M" 1 2 R_0."$N$M"_K"$Z".txt $Z 1.0
		done	
	done
done
#below needs matlab commandline callable and I don't know if -nodisplay works..
./matlab -nodisplay makeAll.m
