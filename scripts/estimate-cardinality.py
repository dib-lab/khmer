#! /usr/bin/env python
"""
Count distinct k-mer in a sequence or a input fa file
,based on hyperloglog algorithm and sha1 hash function. It can have 
small (2-3 %) estimation errors.

% scripts/estimate-cardinality.py <bitsize> <ksize> <filename>

"""

import khmer
import sys
import argparse


default_k=21
default_b=10

def main():
	
	parser = argparse.ArgumentParser()

	parser.add_argument('-k', default=default_k, type=int, help='k-mer size', dest='ksize')
	parser.add_argument('-b', default=default_b, type=int, help='bit size', dest='bsize')
	parser.add_argument('-i', type=str, help='input fasta file', dest='filename')

	args = parser.parse_args()

	k_size = args.ksize
	bit_size= args.bsize
	filename=args.filename

	
	counter=khmer.KmerCardinality(bit_size, k_size)
	counter.consume_fasta(str(filename))
	print 'Estimated k-mers of length',k_size,':',counter.cardinality() 

if __name__ == '__main__':
	main()