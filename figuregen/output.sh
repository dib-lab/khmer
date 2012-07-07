#!/bin/sh

for p in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15
do
   python procresults.py $p/log.txt $p
   Rscript part.R novar.txt var.txt
done
