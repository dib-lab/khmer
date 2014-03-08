#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
"""
simreads.py

Creates short reads from a given genome.
"""
import sys
import random
import math
import screed

def getRandomNucleotide():
   '''
   Returns a nucleotide with even chance for all four.
   '''
   nucleotides = ['A', 'C', 'G', 'T']

   return nucleotides[random.randint(0,3)]

def generateChromosome(length):
   '''
   Generates a simple chromosome with random nucleotides.
   '''
   chrList = []

   for i in range(length):
      chrList.append(getRandomNucleotide())

   return ''.join(chrList)

def generateGenome(genomeSize, chrs):
   '''
   Generates a collection of chromosomes and returns genome as a list of
   chromosome strings.
   '''
   chrSize = int(genomeSize / chrs)
   genome = []
   
   for i in range(chrs):
      genome.append(generateChromosome(chrSize))

   return genome

def writeGenomeFASTA(filename, genome, prefix):
   '''
   Writes the genome to a file in FASTA format.
   '''
   fd = open(filename, "w")

   for i in range(len(genome)):
      fd.write(">" + prefix + ":chr" + str(i+1) + "\n")
      lines = int(math.ceil(len(genome[i])/80.0))
      for j in range(lines):
         fd.write(genome[i][(j*80):((j*80)+80)] + "\n")

def main():

   if len(sys.argv) != 5:
      print "Usage: ./simgenome.py filename genome_size n_chromosomes prefix"
      return

   filename = sys.argv[1]
   genome_size = int(sys.argv[2])
   n_chromosomes = int(sys.argv[3])
   prefix = sys.argv[4]

   genome = generateGenome(genome_size, n_chromosomes)
   writeGenomeFASTA(filename, genome, prefix)

main()
