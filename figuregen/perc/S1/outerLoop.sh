#!/bin/bash
for N in 0 1 2 3 4 5 6 7 8 9; do
	for M in 0 1 2 3 4 5 6 7 8 9; do
		qsub -v localN=$N$M runThis.sh 
	done
done
