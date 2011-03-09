import sys
import khmer
import screed
from screed.fasta import fasta_iter

K=32
HASHTABLE_SIZE = 1e9
N_HT = 4

'''Usage contig.coverage.py <output file> <reference> <reads:>'''

def slidingWindow(sequence, K_size):
    try:
        check = iter(sequence)
    except TypeError:
        raise Exception('ERROR:  Sequence iteration fail.')

    # Total slides needed to consume sequence
    total_windows = (len(sequence) - K_size) + 1

    for i in range(0, total_windows, 1):
        yield sequence[i:i+K_size]


if '__name__==__main__':

    fp = open(sys.argv[1], 'w')
    contig_file = sys.argv[2]
    reads_file = sys.argv[3:]

    ht = khmer.new_counting_hash(K, HASHTABLE_SIZE, N_HT)

    '''consumes contig into hashtable'''
    print '...consuming contig reads from %s' % contig_file
    for n, record in enumerate(fasta_iter(open(contig_file))):
        sequence = record['sequence']
        contig_kmers = slidingWindow(sequence, K)
        for x in contig_kmers:
            if x.find('N') > -1:
                continue
            else:
                ht.consume(x)
    
    '''counts reads into hashtable abundance'''
    for each_file in reads_file:
        read_file = open(each_file, 'r')
        print '...counting in read file from %s' % each_file
        for n1, record1 in enumerate(fasta_iter(read_file)):
            sequence = record1['sequence']
            read_kmers = slidingWindow(sequence, K)
            for kmer in read_kmers:
                if ht.get(kmer) > 0:
                    ht.count(kmer)
        read_file.close()
        
    '''retrieve abundances'''
    d = {}    
    print 'writing output of abundances to %s' % sys.argv[1]
    for n2, record2 in enumerate(fasta_iter(open(contig_file))):
        contig_seq = record2['sequence']
        id = record2['name']
        count_list_forward = []
        contig_kmers = slidingWindow(contig_seq, K)
        for contig_kmer in contig_kmers:
            if contig_kmer.find('N') > -1:
                count_list_forward.append(0)
                continue
            else:
                count_kmer = int(ht.get(contig_kmer)) - 1
                count_list_forward.append(count_kmer)
        for item in count_list_forward:
            print >>fp, '%s' % item

    print 'Hashtable occupancy =', ht.n_occupied()/float(HASHTABLE_SIZE)
