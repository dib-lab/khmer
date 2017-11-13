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
    p.add_argument('--prefix', default='transcriptome')
    args = p.parse_args()

    cg = create_countgraph(args)
    asm = khmer.JunctionCountAssembler(cg)

    tr_fn = '{0}.transcripts.fa'.format(args.prefix)
    orf_fn = '{0}.orfs.fa'.format(args.prefix)
    stats_fn = '{0}.stats.fa'.format(args.prefix)

    with open(tr_fn, 'w') as tr_fp,\
         open(orf_fn, 'w') as orf_fp,\
         open(stats_fn, 'w') as stats_fp:

        kept = 0
        next_contig = 1
        next_orf = 1
        output = set()
        statswriter = csv.DictWriter(stats_fp, delimiter=',',
                                     fieldnames=['read_n', 'action', 'cov',
                                                 'n_junctions', 'contig_n'])

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
                                          'n_junctions': None, 
                                          'contig_n': None})
                elif cov < 30:
                    seq, pos = cg.trim_on_abundance(record.sequence, 3)
                    if len(seq) < args.ksize:
                        continue
                    
                    n_junctions = asm.consume(seq)
                    statswriter.writerow({'read_n': n, 'action': 't', 'cov': cov,
                                          'n_junctions': n_junctions, 
                                          'contig_n': None})
                elif cov == 30:
                    contigs = asm.assemble(record.sequence[:args.ksize])
                    for contig_n, contig in enumerate(contigs):
                        statswriter.writerow({'read_n': n, 'action': 'a', 'cov': cov,
                                              'n_junctions': None, 
                                              'contig_n': (next_contig, contig_n)})
                        tr_fp.write('>contig%d\n%s\n' % (next_contig, contig))
                        next_contig += 1

                        for t in translate(contig):
                            for orf_n, o in enumerate(extract_orfs(t)):
                                if hash(o) not in output:
                                    new = True
                                    output.add(hash(o))
                                    orf_fp.write('>orf%d\n%s\n' % (next_orf, o))
                                    next_orf += 1
                                else:
                                    new = False
                else:
                    statswriter.writerow({'read_n': n, 'action': 's', 'cov': cov,
                                          'n_junctions': None, 
                                          'contig_n': None})

if __name__ == '__main__':
    main()
