# In-progress read-buffering approach to writing out colors to many files
# Basic idea is to buffer some number of reads in memory, then dump them all at once
# Hope that each file acrues, on average, BUFFER_SIZE / NUM_PARTS reads
# ie, if we buffer 1000000 reads, and we have 100000 partitios/colors,
# we should expect the mean buffer size to be 10 reads

import screed
import sys
import argparse
import time

def fastp_iter(filename):
    for record in screed.open(filename, parse_description=False):
        name = record.name
        try:
            name, partition_id = name.rsplit('\t', 1)
        except ValueError:
            print >>sys.stderr, '%%% ERROR: Derp! Is this file partitioned? %%%'
            sys.exit(1)
        # convert name to blast format if necessary
        nname = name.split('|', 2)
        if len(nname) >= 2:
            name = nname[2]
        name = name.split(' ')[0]
        yield name, int(partition_id), record.sequence

class Seq:

    def __init__(self, name, color, seq):
        self.name = name
        self.color = color
        self.seq = seq

    def write(self, fp):
        fp.write('>{}\t{}\n{}\n'.format(self.name, self.color, self.seq))

class ReadBuffer:

    def __init__(self, max_buffers=10000, max_reads=1000000, est_files=100000, output_pref='reads_'):
        self.buffers = {}
        self.buffer_counts = {}
        self.max_buffers = max_buffers
        self.max_reads = max_reads

        self.est_files = est_files
        self.output_pref = output_pref
        self.buffer_flush = self.max_reads / self.est_files

        self.cur_reads = 0
        self.cur_files = 0

    def add_seq(self, seq):
        color = seq.color
        if color in self.buffers:
            count = self.buffer_counts[color]
            self.buffers[color].append(seq)
            self.buffer_counts[color] += 1
            if count > self.buffer_flush:
                self.flush_buffer(color)
                self.del_buffer(color)

        else:
            self.buffers[color] = [seq]
            self.buffer_counts[color] = 1
        self.cur_reads += 1
        if self.cur_reads > self.max_reads:
            self.flush_all()
        if len(self.buffers) > self.max_buffers:
            #self.clean_buffers(2)
            self.flush_all()
    
    def flush_buffer(self, color):
        with open('{}{}.fa'.format(self.output_pref, color), 'a') as outfp:
            for read in self.buffers[color]:
                read.write(outfp)
                self.cur_reads -= 1
            
    def del_buffer(self, color):
        del self.buffer_counts[color]
        del self.buffers[color]

    def flush_all(self):
        print >>sys.stderr, '** reached max buffer size, flushing all to files...'
        for color in self.buffers:
            self.flush_buffer(color)
        colors = self.buffers.keys()
        for color in colors:
            self.del_buffer(color)
        del colors
        assert self.cur_reads == 0

    def clean_buffers(self, cutoff):
        print >>sys.stderr, '** flushing low-abundance buffers...'
        flushed = []
        for color in self.buffers:
            if self.buffer_counts[color] < cutoff:
                self.flush_buffer(color)
                flushed.append(color)
        for color in flushed:
            self.del_buffer(color)
        del flushed

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--buffer_size', dest='buffer_size', type=int)
    parser.add_argument('-e', '--files_estimate', dest='files_estimate', type=int)
    parser.add_argument('-o', '--output_prefix', dest='output_prefix')
    parser.add_argument('-m', '--max_buffers', dest='max_buffers', type=int)
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()

    max_buffers = args.max_buffers
    output_pref = args.output_prefix
    buf_size = args.buffer_size
    est = args.files_estimate
    input_files = args.input_files

    output_buffer = ReadBuffer(max_buffers=max_buffers, max_reads=buf_size, est_files=est, output_pref=output_pref)

    multi_fp = open('{}_multi.fa'.format(output_pref), 'a')
    
    n_reads = 0
    total_t = 0.0
    start_t = time.clock()
    for input_file in args.input_files:
        print >>sys.stderr, '* splitting reads in {}...'.format(input_file)

        current_read = ''
        seen_twice = False

        for name, color, seq in fastp_iter(input_file):
            n_reads += 1
            seq_obj = Seq(name, color, seq)

            if n_reads % 100000 == 0:
                end_t = time.clock()
                batch_t = end_t - start_t
                total_t += batch_t
                print >>sys.stderr, '** processed {} reads from {} [{}s, {}s total]'.format(n_reads, input_file, batch_t, total_t)
                start_t = time.clock()
 
            if name == current_read:
                if not seen_twice:
                    seq_obj.write(multi_fp)
                seen_twice = True
            
            else:
                seen_twice = False
                output_buffer.add_seq(Seq(name,color,seq))
            current_read = name

if __name__ == '__main__':
    main()
