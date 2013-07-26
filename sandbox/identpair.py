import sys
import khmer
from khmer import thread_utils

def filter_fn(a, b, c, d):
    return b, d

n_threads = 8
outfp = open(sys.argv[2], 'w')
tw = thread_utils.PairThreadedWriter(outfp).start()

rparser = khmer.ReadParser(sys.argv[1], n_threads - 1)
threads = thread_utils.start_threads(n_threads - 1,
                                     target=tw.process_fn,
                                     args=(rparser, filter_fn))

tw.join(threads)
