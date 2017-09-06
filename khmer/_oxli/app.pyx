# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from __future__ import print_function
import argparse
import itertools
import json
import os
import sys

from khmer.khmer_args import build_counting_args, create_countgraph
from khmer.khmer_logger import (configure_logging, log_info, log_error,
                                log_warn)

from libcpp cimport bool

from khmer._oxli.graphs cimport Nodegraph, Countgraph

from khmer._oxli.partitioning cimport StreamingPartitioner, Component
from khmer._oxli.partitioning import StreamingPartitioner, Component

from khmer._oxli.parsing cimport BrokenPairedReader, SplitPairedReader, FastxParser, Sequence
from khmer._oxli.parsing import BrokenPairedReader, SplitPairedReader, FastxParser, Sequence
from khmer._oxli.utils cimport _bstring

def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])

cdef class PartitioningApp:

    def __init__(self, args=sys.argv[1:]):
        self.args = self.parse_args(args)
        self.args.write_results = self.args.output_interval > 0

        self.graph = create_countgraph(self.args)
        self.partitioner = StreamingPartitioner(self.graph, tag_density=self.args.tag_density)

    def parse_args(self, args):
        parser = build_counting_args(descr='Partition a sample',
                                     citations=['counting', 'SeqAn'])
        parser.add_argument('--output-dir', default='partitioned')
        parser.add_argument('samples', nargs='+')
        parser.add_argument('--save', action='store_true', default=False)
        parser.add_argument('--pairing-mode', 
                            choices=['split', 'interleaved', 'single'],
                            default='split')
        parser.add_argument('-Z', dest='norm', default=10, type=int)
        parser.add_argument('--output-interval', default=0, type=int)
        parser.add_argument('--tag-density', default=None, type=int)
        
        return parser.parse_args(args)

    def write_results(self, folder, n, new_kmers):
        filename = os.path.join(folder, '{0}.csv'.format(n))
        print('# {0}: {1} tags, {2} components.'.format(n, self.partitioner.n_tags, 
                                                        self.partitioner.n_components))
        print('  writing results to file -> {0}'.format(filename))
        self.partitioner.write_components(filename)
        with open(os.path.join(folder, 'global.csv'), 'a') as fp:
            fp.write('{0}, {1}, {2}, {3}\n'.format(n, self.partitioner.n_components,
                                                 self.partitioner.n_tags, new_kmers))
        cov_filename = os.path.join(folder, '{0}.coverage.csv'.format(n))
        self.partitioner.write_component_coverage(cov_filename)

    def prep_results_dir(self):
        try:
            os.mkdir(self.args.output_dir)
        except OSError as e:
            pass

        if self.args.save:
            self.args.save = os.path.join(self.args.output_dir, 'partitioner')

    def write_meta(self, n_sequences, total_kmers):
        meta = {'samples': self.args.samples,
                'pairing': self.args.pairing_mode,
                'K': self.args.ksize,
                'tag-density': self.partitioner.tag_density,
                'n_sequences': n_sequences,
                'n_unique_kmers': total_kmers}
        if self.args.save:
            meta['partitioner'] = self.args.save

        with open(os.path.join(self.args.output_dir, 'meta'), 'w') as fp:
            json.dump(meta, fp, indent=4)

    def run(self):

        self.prep_results_dir()

        if self.args.pairing_mode == 'split':
            samples = list(grouper(2, self.args.samples))
            for pair in samples:
                if len(pair) != 2:
                    raise ValueError('Must have even number of samples!')
        else:
            samples = self.args.samples
        
        cdef int n
        cdef int n_sequences = 0
        cdef bool paired
        cdef Sequence first, second
        cdef int new_kmers = 0
        cdef int total_kmers = 0
        cdef int print_interval = self.args.output_interval if self.args.write_results else 10000
        last = 0
        for group in samples:
            if self.args.pairing_mode == 'split':
                sample_name = '{0}.{1}'.format(group[0], group[1])
                print('== Starting ({0}) =='.format(sample_name))
                reader = SplitPairedReader(FastxParser(group[0]),
                                           FastxParser(group[1]),
                                           min_length=self.args.ksize)
            else:
                sample_name = group
                print('== Starting {0} =='.format(sample_name))
                reader = BrokenPairedReader(FastxParser(group), min_length=self.args.ksize)
            for n, paired, first, second in reader:

                if n % print_interval == 0:
                    print (n, self.partitioner.n_components, self.partitioner.n_tags)
                if self.args.write_results and n > 0 and n % self.args.output_interval == 0:
                    self.write_results(self.args.output_dir, last+n, new_kmers)
                    total_kmers += new_kmers
                    new_kmers = 0
                if paired:
                    new_kmers += self.partitioner.consume_pair(first.sequence,
                                                               second.sequence)
                else:
                    new_kmers += self.partitioner.consume(first.sequence)
            last = n
            n_sequences += last
            if self.args.write_results:
                self.write_results(self.args.output_dir, last, new_kmers)
                total_kmers += new_kmers
                new_kmers = 0

        if self.args.save:
            self.partitioner.save(self.args.save)

        self.write_meta(n_sequences, total_kmers)

        return self.partitioner


cdef class DynamicPartitioning(PartitioningApp):

    def run(self):
        pass
