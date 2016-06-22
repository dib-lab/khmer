simple-genome-reads.fa
----------------------

* generated from https://github.com/dib-lab/nullgraph,
  commit 2146df2017d0dabc22a5a

* commands::

   make-random-genome.py -l 1000 -s 1 > simple-genome.fa
   make-reads.py -S 1 -e .01 -r 100 -C 100 simple-genome.fa > simple-genome-reads.fa
