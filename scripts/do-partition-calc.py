#! /usr/bin/env python
import khmer, sys
from optparse import OptionParser

### defaults

CHECKPOINT_PERIOD=1000000               # 1m
K=32
HASHTABLE_SIZE=15

###

parser = OptionParser()
parser.add_option("-k", "--k", dest="ksize", help="k-mer size", default=K,
                  type=int)
parser.add_option("-s", "--size", dest="htsize4", help="hashtable size",
                  default=HASHTABLE_SIZE, type=int)
parser.add_option("-c", "--checkpoint", dest="checkpoint",
                  help="checkpoint period", 
                  default=CHECKPOINT_PERIOD, type=int)
parser.add_option("--load", dest="load_partitionmap_from",
              help="load partitionmap data from <file>",
                  default=None)

###

def make_reporting_fn(ht, filename, period):

    def _report(name, count1, count2,
                ht=ht, f=filename, checkpoint_period=period):
        print name, count1, count2
        if name == 'do_truncated_partition/read' and \
               count1 % checkpoint_period == 0:
            ht.save_partitionmap(f)

    return _report

def main():
    (options, args) = parser.parse_args()

    (infile,) = args
    outfile = infile + '.part'

    k = options.ksize
    hashtable_size = 4**options.htsize4 + 1
    checkpoint_period = options.checkpoint

    if k > 32:
        raise Exception, "invalid k; must be <= 32"

    ###

    print 'making hashtable: k=%d, hashtable size=%.2f bn' % (k,
                                                 hashtable_size / float(1e9))
    ht = khmer.new_hashtable(k, hashtable_size)

    ###

    if options.load_partitionmap_from:
        print '** loading partitionmap data from', \
              options.load_partitionmap_from
        f = options.load_partitionmap_from
        ht.load_partitionmap(f)

    ###

    report_fn = make_reporting_fn(ht, outfile, checkpoint_period)
    
    n_partitions = ht.do_truncated_partition(infile, outfile, report_fn)
    print n_partitions, 'partitions kept'

    ht.save_partitionmap(outfile + '.end.pmap')

if __name__ == '__main__':
    main()
