#! /usr/bin/env python
import khmer, sys

CHECKPOINT_PERIOD=1000000
K=32
HASHTABLE_SIZE=int(4**15)+1

def make_reporting_fn(ht, filename, period):

    def _report(name, count1, count2,
                ht=ht, f=filename, checkpoint_period=period):
        print name, count1, count2
        if name == 'do_truncated_partition/read' and \
               count1 % checkpoint_period == 0:
            ht.save_checkpoint(f + '.pmap',
                               f + '.surrender')

    return _report

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    k = K
    hashtable_size = HASHTABLE_SIZE
    checkpoint_period = CHECKPOINT_PERIOD

    print 'making hashtable: k=%d, hashtable size=%dbn' % (k,
                                                          hashtable_size / 1e9)
    ht = khmer.new_hashtable(k, hashtable_size)

    report_fn = make_reporting_fn(ht, outfile, checkpoint_period)
    
    n_partitions = ht.do_truncated_partition(infile, outfile, report_fn)
    print n_partitions, 'partitions kept'

    ht.save_checkpoint(infile + '.pmap.end',
                       outfile + '.surrender.end')

if __name__ == '__main__':
    main()

    
