#!/bin/sh -login
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=7gb
#PBS -j oe

cd /mnt/scratch/pelljaso/perc/S1

for Z in 4 5 6 7 8 9 10 11 12; do
  ./ppt 0."$localN" 1 100 R_0."$localN"_K"$Z".txt $Z 1.0
done
