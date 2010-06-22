import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fastq import fastq_iter
import glob
import khmer
import random

def getRandomNucleotide():
   '''
   Returns a nucleotide with even chance for all four.
   '''
   nucleotides = ['A', 'C', 'G', 'T']

   return nucleotides[random.randint(0,3)]


def generateKmer(length):
   '''
   Generates a simple chromosome with random nucleotides.
   '''
   chrList = []

   for i in range(length):
      chrList.append(getRandomNucleotide())

   return ''.join(chrList)

def complement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
   complseq = [complement[base] for base in seq]
   return complseq

def reverse_complement(seq):
   seq = list(seq)
   seq.reverse()
   return ''.join(complement(seq))

K = 17
THRESHOLD = 3
MINLENGTH = 50

hasher = khmer.HashtableIntersect(K, *khmer.PRIMES_8b)

files = glob.glob(sys.argv[1])

for filename in files:
    for n, record in enumerate(fastq_iter(open(filename))):
        if n % 10000 == 0:
            print>>sys.stderr, '...', n
        if n > 10**6*4:
            break
        
        sequence = record['sequence']
        if 'N' in sequence:
            continue

        revComp = reverse_complement(sequence)
        if revComp < sequence:
            hasher.consume(revComp)
        else:
            hasher.consume(sequence)

falseCount = 0
numTried = 1000

for i in range(numTried):
   kmer = generateKmer(K)
   if hasher.get_min_count(kmer) > 0:
      falseCount += 1

print falseCount
