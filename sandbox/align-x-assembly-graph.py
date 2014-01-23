#! /usr/bin/python
'''

necessary components:
- label to sequence dictionary
- simple smith-waterman sequence aligner

basic steps:
1) generate screed database
2) consume and tag database
3) store label to sequence dict
4) sweep labels from queries (in parallel ideally)
5) align queries to label->database sequences

'''

import screed
from screed import ScreedDB
import khmer
from khmer import LabelHash

from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2

import argparse
import sys
import os

DEFAULT_K = 21
DEFAULT_HTSIZE=1e8
DEFAULT_N_HT=4

TRAVERSE=1

ALN_MAT=matlist.blosum62
GAP_OPEN=-7
GAP_EXTEND=-0.5

def index_db(dbfile):
    screed_db_name = '{}_screed'.format(dbfile) 
    if os.path.exists(screed_db_name):
        print >>sys.stderr, 'data file {} already exists, skipping indexing...'.format(screed_db_name)
    else:
        print >>sys.stderr, 'indexing {} using screed...'.format(dbfile)
        screed.read_fasta_sequences(dbfile)
    return screed_db_name

def align_sequences(qseq, dseq):
    alignments = pairwise2.align.localdx(qseq, dseq, ALN_MAT)
    if alignments:
        top = alignments[0]
        aln_qseq, aln_dseq, score, begin, end = top
        return (score, begin, end)
    else:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--ksize', dest='ksize', type=int, default=DEFAULT_K)
    parser.add_argument('-x', '--htize', dest='htsize', type=float, default=DEFAULT_HTSIZE)
    parser.add_argument('-q', '--query', dest='query')
    parser.add_argument('-d', '--db', dest='database')
    parser.add_argument('-o', '--output_name', dest='output_name')
    args = parser.parse_args()

    # should probably check for file existnece as well
    if not args.query:
        print >>sys.stderr, 'No query sequence specified, exiting...'
        sys.exit()
    if not args.database:
        print >>sys.stderr, 'No database sequence specified, exiting...'
        sys.exit()

    try:
        outfp = open(args.output_name, 'wb')
    except IOError as e:
        print >>sys.stderr, 'error opening output file'
        print e
        sys.exit()

    outfp.write('qname,sname,slabel,minident,ntags\n')

    screed_db_name = index_db(args.database)
    db = ScreedDB(screed_db_name)
    
    lh = LabelHash(args.ksize, args.htsize, DEFAULT_N_HT)
    
    print >>sys.stderr, 'consuming and labeling database...'
    for n, key in enumerate(db.keys()):
        record = db[key]
        seq = str(record.sequence)
        
        if n % 25000 == 0:
            print >>sys.stderr, '...consumed and labeled {} sequences'.format(n)
            print >>sys.stderr, '......label: {}'.format(n)
            #print seq, sid, type(seq), type(sid)

        lh.consume_sequence_and_tag_with_labels(seq, n)

    print >>sys.stderr, 'beginning query sequence alignment...'
    for n, record in enumerate(screed.open(args.query)):
        if n % 25000 == 0:
            print >>sys.stderr, '...aligned {} query sequences'.format(n)
        
        name = record.name
        seq = record.sequence
       
        tags = lh.sweep_tag_neighborhood(seq, args.ksize)

        labels_by_tag = {}
        tags_by_label = {}
        labels = set()
        for tag in tags:
            tmp = lh.get_tag_labels(tag)
            labels.update(tmp) 
            labels_by_tag[tag] = tmp
            for label in tmp:
                if label not in tags_by_label:
                    tags_by_label[label] = set()
                tags_by_label[label].add(tag)
        
        for label in tags_by_label:
            tags = tags_by_label[label]
            dbrecord = db.loadRecordByIndex(label)
            dbseq_name = str(dbrecord.name)
            mident = float(args.ksize * len(tags)) / len(seq)
            nkmers = len(tags)
            outfp.write('{q},{d},{id},{s},{k}\n'.format(q=name, d=dbseq_name, id=label,
                                                        s=mident, k=nkmers))
             
        ''' 
        labels = lh.sweep_label_neighborhood(seq, TRAVERSE)
        if labels:
            for label in labels:
                dbrecord = db.loadRecordByIndex(label)
                dbseq = str(dbrecord.sequence)
                aln = align_sequences(seq, dbseq)
                if aln:
                    score, begin, end = aln
                    length = max(end-begin, begin-end)

                    dbname = str(dbrecord.name)
                    outfp.write('{q},{d},{id},{s},{l},{b},{e}\n'.format(q=name, d=dbname, id=label,
                                                                      s=score, l=length, b=begin, e=end))
        '''

if __name__ == '__main__':
    main() 
