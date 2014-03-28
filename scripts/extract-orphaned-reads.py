import screed
import os
import sys
import argparse
from collections import deque

def validpair(r0, r1):
    return r0.name[-1] == "1" and \
        r1.name[-1] == "2" and \
        r0.name[0:-1] == r1.name[0:-1]

def write(record, fp):
    if hasattr(record, 'accuracy'):
        fp.write('@{name}\n{seq}\n+\n{acc}\n'.format(name=record.name,
                                seq=record.sequence,
                                acc=record.accuracy))
    else:
        fp.write('>{name}\n{seq}\n'.format(name=record.name,
                                      seq=record.sequence))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_filenames', nargs='+')
    args = parser.parse_args()

    Q = deque()

    for filename in args.input_filenames:

        pe_fn = os.path.basename(filename) + '.pe'
        se_fn = os.path.basename(filename) + '.se'
        with open(pe_fn, 'wb') as pe_fp, open(se_fn, 'wb') as se_fp:
            n_se = 0
            n_pe = 0
            for n, record in enumerate(screed.open(filename)):
                if n % 100000 == 0:
                    print >>sys.stderr, n, 'processed...'

                Q.append(record)
                if len(Q) == 2:
                    
                    if validpair(Q[0], Q[-1]):
                        n_pe += 2
                        write(Q.popleft(), pe_fp)
                        write(Q.popleft(), pe_fp)
                    else:
                        n_se += 1
                        write(Q.popleft(), se_fp)

        assert n == n_pe+n_se-1
        print filename, '-- {npe} paired, {nse} orphaned'.format(
                                                        npe=n_pe,
                                                        nse=n_se)

if __name__ == '__main__':
    main()
