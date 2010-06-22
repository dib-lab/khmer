#!/bin/sh

PYTHONPATH=../python python fasta-filter.py /scratch/jason/1.fa 8000000011 > /scratch/jason/8000000011.fa

j=8000000011

for i in 8000000111 8000000257 8000000321 8000000381
do
   PYTHONPATH=../python python fasta-filter.py /scratch/jason/$j.fa $i > /scratch/jason/$i.fa

   j=$i
done
