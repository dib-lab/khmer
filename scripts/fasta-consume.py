import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fasta import fasta_iter
import glob
import khmer

def complement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
   complseq = [complement[base] for base in seq]
   return complseq

def reverse_complement(seq):
   seq = list(seq)
   seq.reverse()
   return ''.join(complement(seq))


def getKmer(seq):
    revComp = reverse_complement(seq)
    if revComp < sequence:
        return revComp
    else:
        return seq


K = 17
THRESHOLD = 3
MINLENGTH = 50

#hasher = khmer.HashtableIntersect(K, *khmer.PRIMES_4b)
hasher = khmer.new_hashtable(K, khmer.PRIMES_8b[0])

files = glob.glob(sys.argv[1])

for filename in files:
   hasher.consume_fasta(filename)

'''
for filename in files:
    for n, record in enumerate(fasta_iter(open(filename))):
        if n % 10000 == 0:
            print>>sys.stderr, '...', n
        #if n > 10**5:
        #    break
        
        sequence = record['sequence']
        if 'N' in sequence:
            continue
        
        hasher.consume(getKmer(sequence))
'''

for filename in files:
    for n, record in enumerate(fasta_iter(open(filename))):
        if n % 10000 == 0:
            print>>sys.stderr, '2...', n
        #if n > 10**5:
        #    break
        
        sequence = record['sequence']
        if 'N' in sequence:
            continue

        readAbund = [0]*98
        for i in range(len(sequence) - K + 1):
            kmer = getKmer(sequence[i:(i+K)])
            readAbund[i] += hasher.get_min_count(kmer)
        
        start = 0
        for i in range(len(readAbund)):
            if readAbund[i] >= THRESHOLD:
                break
            else:
                start += 1
            
        stop = len(readAbund)-1    
        for i in range(len(readAbund)-1, -1, -1):
            if readAbund[i] >= THRESHOLD:
                break
            else:
                stop -= 1

        if start == len(readAbund) or (stop - start + K) < MINLENGTH:
            continue
        else:
            mySeq = sequence[start:stop+K]
            print ">" + record['name']
            print mySeq
            #print start, stop
            #print readAbund
