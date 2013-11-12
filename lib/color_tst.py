import khmer
import screed

def reverse_comp(s):
    ret = ''
    for i in range(len(s)-1,-1,-1):
        c = s[i]
        if c == 'A':
            ret += 'T'
        elif c == 'T':
            ret += 'A'
        elif c == 'G':
            ret += 'C'
        else:
            ret += 'G'
    return ret

ht = khmer.new_hashbits(20,1e8,4)
print '#' * 200
ht.consume_fasta_and_tag_with_colors('../tests/test-data/test-reads.fa')
#print ht.sweep_sequence_for_colors('CACACACGGACATCGGAGAGAGGCTGAGACAGCGAGACACACAGAGACAGAGCGGAGAGGGCACAGACAGACAAGAGCATGAGAGATCGGCAGAGCGGTG', False, False)
#print ht.sweep_sequence_for_colors('CGCCGTAGTCGTACTGGTTCTCCTCCGTGTACTCGTGCGCTGCCTCCACCTCTGGGCTGCTCATGCCCTCCATGTGACCTTCAGGCATGCCCTCGGAGAT', False, False)
#print ht.sweep_sequence_for_colors('GGAGAGCCTGGGGCCAAGCCCGAGGGCATGCCTGAAGGTCACATGGAGGGCATGAGCAGCCCAG', False, False)
#print ht.sweep_sequence_for_colors('TTTTTTGAATACGTTTAGTTAATATTTGTACTTCAATTAATAAAAATTTGCTATAATTTTTCCATTATCGCCAGTCACTCGCGTGATATAGGAAAAGGTT', False, False)
#print ht.sweep_sequence_for_colors('AAGCAGTGGTATCAACGCAGAGTACGCGGGGACTCTGTCGCTGCTCCTCTAGCACAGAGAGCCAGAGACGGCTTACAGCAGCAGCATCATATAGCCTC', False, False)

t0 = 'CCATGTAGCGCCGCACACCTTTGTAGGTGTTGTAATAATCTTCGATGACTTTCTTCGCTTCCTGACGGCTTATGCC'
t1 = 'ACCGCGCGCGAATCGACGGTTGTCAGCCAAAGGCGTTCAACACCAGCACCGCCCTTAAGCCGCCCGCCCGCCGCCC'
N = 1000

for n, record in enumerate(screed.open('../tests/test-data/test-reads.fa')):
    if n > N:
        break
    print '*' * 40
    seq = record.sequence
    print seq
    colors = ht.sweep_sequence_for_colors(seq, False, False)
    print 'colors from sweep:', colors
    tags = ht.get_all_tags(seq)
    print 'tags from get_all_tags:', tags
    print 'colors from get_tag_colors:'
    t_colors = set()
    for tag in tags:
        t_colors.update(ht.get_tag_colors(tag))
    print t_colors
    assert len(t_colors) == len(colors)

'''
file_pointers = {}
for n, record in enumerate(screed.open('/w/2013-lamprey/syn_part/syn.sweep.fa')):
    if n >= N:
        break
    if n % 1000 == 0:
        print '...processed {} reads'.format(n)
    colors = ht.sweep_sequence_for_colors(record.sequence, False, False)
    for c in colors:
        if c in file_pointers.viewkeys():
            file_pointers[c].write('>{}\n{}\n'.format(record.name, record.sequence))
        else:
            file_pointers[c] = open('color_{}.fa'.format(c), 'wb')
            file_pointers[c].write('>{}\n{}\n'.format(record.name, record.sequence))\
'''
'''
ht = khmer.new_hashbits(25, 1e9,4)
ht.consume_partitioned_fasta_and_tag_with_colors('/w/2013-lamprey/test.fp')

for n, record in enumerate(screed.open('/w/lamprey-mrnaseq/reads/single/L82-a.fq.gz')):
    if n >= N:
        break
    colors = ht.sweep_sequence_for_colors(record.sequence, False,  False)
    if colors:
        print colors
'''
