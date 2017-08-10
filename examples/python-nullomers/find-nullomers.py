#! /usr/bin/env python
import khmer
import argparse
import screed

K = 11
SIZE = 4**K                               # important: use exact counting.

def main():
    p = argparse.ArgumentParser()
    p.add_argument('contigs', nargs='+')
    args = p.parse_args()

    assert K % 2 == 1, "K must be odd"
    
    print('allocating lots of memory for exact counts: {} bytes'.format(8*SIZE))
    # use Countgraph because we want to use the predictable hash function
    # that only goes up to 4**K.
    ct = khmer.Countgraph(K, SIZE, 1)

    for filename in args.contigs:
        print('consuming {}'.format(filename))
        ct.consume_seqfile(filename)
    print('...done!')

    print('Iterating over all {}-mers'.format(K))

    # for large K, this is going to end up producing a massive amount of
    # output...
    for i in range(SIZE):
        print(ct.reverse_hash(i), ct.get(i))


if __name__ == '__main__':
    main()
