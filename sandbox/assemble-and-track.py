#! /usr/bin/env python
import csv
import screed
import khmer
from khmer.khmer_args import build_counting_args, create_countgraph
import argparse
import sys


DEFAULT_COV = 20
THRESH2 = 30


dna_to_aa={'TTT':'F','TTC':'F', 'TTA':'L','TTG':'L',
                'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y', 'TAA':'*','TAG':'*','TGA':'*',
                'TGT':'C','TGC':'C', 'TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H', 'CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I', 'ATG':'M',
                'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N', 'AAA':'K','AAG':'K',
                'AGT':'S','AGC':'S', 'AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D', 'GAA':'E','GAG':'E',
                'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}


__complementTranslation = { "A": "T", "C": "G", "G": "C", "T": "A", "N": "N" }
def complement(s):
    """
    Return complement of 's'.
    """
    c = "".join(__complementTranslation[n] for n in s)
    return c


def reverse(s):
    """
    Return reverse of 's'.
    """
    r = "".join(reversed(s))
    return r


def peptides(seq, start):
    for i in range(start, len(seq), 3):
        yield dna_to_aa.get(seq[i:i+3], "X")


def translate(seq):
    for i in range(3):
        pep = peptides(seq, i)
        yield "".join(pep)

    revcomp = reverse(complement((seq)))
    for i in range(3):
        pep = peptides(revcomp, i)
        yield "".join(pep)


def extract_orfs(pepseq, min_length=99):
    peplist = [ x for x in pepseq.split('*') if len(x) >= min_length ]
    for x in peplist:
        yield x


def main():
    p = build_counting_args(descr='Streaming assembly with tracking info')
    p.add_argument('fastq_files', nargs='+')
    p.add_argument('-o', type=argparse.FileType('w'),
                   default='assembly-stats.csv')
    args = p.parse_args()

    cg = create_countgraph(args)

    kept = 0
    hdn = khmer.HashSet(args.ksize)
    lh = khmer._GraphLabels(cg)
    next_label = 1
    next_orf = 1
    output = set()
    statswriter = csv.DictWriter(args.o, delimiter=',',
                                 fieldnames=['read_n', 'action', 'cov', 'n_hdn',
                                             'contig_n', 'orf_n', 'new'])

    for filename in args.fastq_files:
        for n, record in enumerate(screed.open(filename)):
            if n and n % 10000 == 0:
                print('...', n, file=sys.stderr)

            if len(record.sequence) < args.ksize:
                continue

            cov, _, _ = cg.get_median_count(record.sequence)
            if cov < 20:
                kept += 1
                cg.consume(record.sequence)
                statswriter.writerow({'read_n': n, 'action': 'c', 'cov': cov,
                                      'n_hdn': None, 'contig_n': None, 
                                      'orf_n': None, 'new': None})
            elif cov < 30:
                #print('intermediate', next_label, file=sys.stderr)
                seq, pos = cg.trim_on_abundance(record.sequence, 3)
                if len(seq) < args.ksize:
                    continue
                
                cg.consume(seq)
                hdn = cg.find_high_degree_nodes(seq)
                lh.label_across_high_degree_nodes(seq, hdn, next_label)
                next_label += 1
                statswriter.writerow({'read_n': n, 'action': 'l', 'cov': cov,
                                      'n_hdn': len(hdn), 'contig_n': None, 
                                      'orf_n': None, 'new': None})
            elif cov == 30:
                contigs = lh.assemble_labeled_path(record.sequence[:args.ksize])
                for contig_n, contig in enumerate(contigs):
                    statswriter.writerow({'read_n': n, 'action': 'a', 'cov': cov,
                                          'n_hdn': None, 'contig_n': contig_n, 
                                          'orf_n': None, 'new': None})
                    for t in translate(contig):
                        for orf_n, o in enumerate(extract_orfs(t)):
                            if hash(o) not in output:
                                new = True
                                output.add(hash(o))
                                print('>orf%d\n%s' % (next_orf, o))
                                next_orf += 1
                            else:
                                new = False
                            statswriter.writerow({'read_n': n, 'action': 'a', 'cov': cov,
                                                  'n_hdn': None, 'contig_n': contig_n, 
                                                  'orf_n': orf_n, 'new': new})

if __name__ == '__main__':
    main()
