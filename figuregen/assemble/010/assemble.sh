#!/bin/sh

for i in msb2.group*.fa
do
   ABYSS -k33 $i -o contigs.$i
done
