#! /usr/bin/env python
import khmer
import argparse
import screed

K = 11
SIZE = 4**K

def main():
    p = argparse.ArgumentParser()
    p.add_argument('contigs', nargs='+')
    args = p.parse_args()

    assert K % 2 == 1, "K must be odd"
    
    print('allocating lots of memory: {}'.format(SIZE))
    ct = khmer.Countgraph(K, SIZE, 1)

    for filename in args.contigs:
        print('consuming {}'.format(filename))
        ct.consume_seqfile(filename)
    print('...done!')

    print('Iterating over all {}-mers'.format(K))

    for i in range(SIZE):
        print(ct.reverse_hash(i), ct.get(i))


if __name__ == '__main__':
    main()
