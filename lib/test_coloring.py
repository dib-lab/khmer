import khmer
import screed


ht = khmer.new_hashbits(20,1e8,4)
print '#' * 200
ht.consume_fasta_and_tag_with_colors('/w/2013-lamprey/syn_part/syn.trinity.fasta')
#print ht.sweep_sequence_for_colors('CACACACGGACATCGGAGAGAGGCTGAGACAGCGAGACACACAGAGACAGAGCGGAGAGGGCACAGACAGACAAGAGCATGAGAGATCGGCAGAGCGGTG', False, False)
#print ht.sweep_sequence_for_colors('CGCCGTAGTCGTACTGGTTCTCCTCCGTGTACTCGTGCGCTGCCTCCACCTCTGGGCTGCTCATGCCCTCCATGTGACCTTCAGGCATGCCCTCGGAGAT', False, False)
#print ht.sweep_sequence_for_colors('GGAGAGCCTGGGGCCAAGCCCGAGGGCATGCCTGAAGGTCACATGGAGGGCATGAGCAGCCCAG', False, False)
#print ht.sweep_sequence_for_colors('TTTTTTGAATACGTTTAGTTAATATTTGTACTTCAATTAATAAAAATTTGCTATAATTTTTCCATTATCGCCAGTCACTCGCGTGATATAGGAAAAGGTT', False, False)
#print ht.sweep_sequence_for_colors('AAGCAGTGGTATCAACGCAGAGTACGCGGGGACTCTGTCGCTGCTCCTCTAGCACAGAGAGCCAGAGACGGCTTACAGCAGCAGCATCATATAGCCTC', False, False)

N=1000000000

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

ht = khmer.new_hashbits(25, 1e9,4)
ht.consume_partitioned_fasta_and_tag_with_colors('/w/2013-lamprey/test.fp')

for n, record in enumerate(screed.open('/w/lamprey-mrnaseq/reads/single/L82-a.fq.gz')):
    if n >= N:
        break
    colors = ht.sweep_sequence_for_colors(record.sequence, False,  False)
    if colors:
        print colors

