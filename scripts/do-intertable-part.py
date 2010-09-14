import sys
import khmer

infile = sys.argv[1]
khmer.do_intersection_partition(19, 250000013, #1000000007,
                                infile, infile + '.ipart')
