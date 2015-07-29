#
# vim: set encoding=utf-8
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

from __future__ import unicode_literals

import sys
import os
import argparse
from argparse import _VersionAction

import screed
import khmer
from khmer import extract_countinghash_info, extract_nodegraph_info
from khmer import __version__
from .utils import print_error
from .khmer_logger import log_info


DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = 1e6
DEFAULT_N_THREADS = 1


class _VersionStdErrAction(_VersionAction):

    def __call__(self, parser, namespace, values, option_string=None):
        version = self.version
        if version is None:
            version = parser.version
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stderr)
        parser.exit()


class ComboFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    pass


def build_graph_args(descr=None, epilog=None, parser=None):
    """Build an ArgumentParser with args for bloom filter based scripts."""
    if parser is None:
        parser = argparse.ArgumentParser(description=descr, epilog=epilog,
                                         formatter_class=ComboFormatter)

    parser.add_argument('--version', action=_VersionStdErrAction,
                        version='khmer {v}'.format(v=__version__))

    parser.add_argument('--ksize', '-k', type=int, default=DEFAULT_K,
                        help='k-mer size to use')

    parser.add_argument('--n_tables', '-N', type=int,
                        default=DEFAULT_N_TABLES,
                        help='number of k-mer counting tables to use')
    parser.add_argument('-U', '--unique-kmers', type=int, default=0,
                        help='approximate number of unique kmers in the input'
                             ' set')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--max-tablesize', '-x', type=float,
                       default=DEFAULT_MAX_TABLESIZE,
                       help='upper bound on tablesize to use; overrides ' +
                       '--max-memory-usage/-M.')
    group.add_argument('-M', '--max-memory-usage', type=float,
                       help='maximum amount of memory to use for data ' +
                       'structure.')

    return parser


def build_counting_args(descr=None, epilog=None):
    """Build an ArgumentParser with args for countinggraph based scripts."""
    parser = build_graph_args(descr=descr, epilog=epilog)
    parser.hashtype = 'countgraph'

    return parser


def build_nodegraph_args(descr=None, epilog=None, parser=None):
    """Build an ArgumentParser with args for nodegraph based scripts."""
    parser = build_graph_args(descr=descr, epilog=epilog, parser=parser)
    parser.hashtype = 'nodegraph'

    return parser

# add an argument for loadhash with warning about parameters


def add_loadgraph_args(parser):

    class LoadAction(argparse.Action):

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, values)

            if getattr(namespace, 'ksize') != DEFAULT_K or \
               getattr(namespace, 'n_tables') != DEFAULT_N_TABLES or \
               getattr(namespace, 'max_tablesize') != DEFAULT_MAX_TABLESIZE:
                if values:
                    print_error('''
** WARNING: You are loading a saved k-mer table from
** {hashfile}, but have set k-mer table parameters.
** Your values for ksize, n_tables, and tablesize
** will be ignored.'''.format(hashfile=values))

            if hasattr(parser, 'hashtype'):
                info = None
                if parser.hashtype == 'nodegraph':
                    info = extract_nodegraph_info(
                        getattr(namespace, self.dest))
                elif parser.hashtype == 'countgraph':
                    info = extract_countinghash_info(
                        getattr(namespace, self.dest))
                if info:
                    K = info[0]
                    x = info[1]
                    n = info[2]
                    setattr(namespace, 'ksize', K)
                    setattr(namespace, 'n_tables', n)
                    setattr(namespace, 'max_tablesize', x)

    parser.add_argument('-l', '--loadtable', metavar="filename", default=None,
                        help='load a precomputed k-mer table from disk',
                        action=LoadAction)


def calculate_tablesize(args, hashtype, multiplier=1.0):
    if hashtype not in ('countgraph', 'nodegraph'):
        raise ValueError("unknown graph type: %s" % (hashtype,))

    if args.max_memory_usage:
        if hashtype == 'countgraph':
            tablesize = args.max_memory_usage / args.n_tables / \
                float(multiplier)
        elif hashtype == 'nodegraph':
            tablesize = 8. * args.max_memory_usage / args.n_tables / \
                float(multiplier)
    else:
        tablesize = args.max_tablesize

    return tablesize


def create_nodegraph(args, ksize=None, multiplier=1.0):
    if ksize is None:
        ksize = args.ksize
    if ksize > 32:
        print_error("\n** ERROR: khmer only supports k-mer sizes <= 32.\n")
        sys.exit(1)

    tablesize = calculate_tablesize(args, 'nodegraph', multiplier)
    return khmer.Hashbits(ksize, tablesize, args.n_tables)


def create_countgraph(args, ksize=None, multiplier=1.0):
    if ksize is None:
        ksize = args.ksize
    if ksize > 32:
        print_error("\n** ERROR: khmer only supports k-mer sizes <= 32.\n")
        sys.exit(1)

    tablesize = calculate_tablesize(args, 'countgraph', multiplier=multiplier)
    return khmer.CountingHash(ksize, tablesize, args.n_tables)


def report_on_config(args, hashtype='countgraph'):
    """Print out configuration.

    Summarize the configuration produced by the command-line arguments
    made available by this module.
    """
    from khmer.utils import print_error
    if hashtype not in ('countgraph', 'nodegraph'):
        raise ValueError("unknown graph type: %s" % (hashtype,))

    tablesize = calculate_tablesize(args, hashtype)

    print_error("\nPARAMETERS:")
    print_error(" - kmer size =    {0} \t\t(-k)".format(args.ksize))
    print_error(" - n tables =     {0} \t\t(-N)".format(args.n_tables))
    print_error(
        " - max tablesize = {0:5.2g} \t(-x)".format(tablesize)
    )
    print_error("")
    if hashtype == 'countgraph':
        print_error(
            "Estimated memory usage is {0:.2g} bytes "
            "(n_tables x max_tablesize)".format(
                args.n_tables * tablesize))
    elif hashtype == 'nodegraph':
        print_error(
            "Estimated memory usage is {0:.2g} bytes "
            "(n_tables x max_tablesize / 8)".format(args.n_tables *
                                                    tablesize / 8)
        )

    print_error("-" * 8)

    if DEFAULT_MAX_TABLESIZE == tablesize and \
       not getattr(args, 'loadtable', None):
        print_error('''\

** WARNING: tablesize is default!
** You probably want to increase this with -M/--max-memory-usage!
** Please read the docs!
''')


def add_threading_args(parser):
    """Add option for threading to options parser."""
    parser.add_argument('--threads', '-T', default=DEFAULT_N_THREADS, type=int,
                        help='Number of simultaneous threads to execute')

_algorithms = {
    'software': 'MR Crusoe et al., '
    '2014. http://dx.doi.org/10.6084/m9.figshare.979190',
    'diginorm': 'CT Brown et al., arXiv:1203.4802 [q-bio.GN]',
    'streaming': 'Q Zhang, S Awad, CT Brown, '
    'https://dx.doi.org/10.7287/peerj.preprints.890v1',
    'graph': 'J Pell et al., http://dx.doi.org/10.1073/pnas.1121464109',
    'counting': 'Q Zhang et al., '
    'http://dx.doi.org/10.1371/journal.pone.0101271',
    'sweep': 'C Scott, MR Crusoe, and CT Brown, unpublished',
    'SeqAn': 'A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11',
    'hll': 'Irber and Brown, unpublished'
}


def info(scriptname, algorithm_list=None):
    """Print version and project info to stderr."""
    import khmer

    log_info("\n|| This is the script {name} in khmer.\n"
             "|| You are running khmer version {version}",
             name=scriptname, version=khmer.__version__)
    log_info("|| You are also using screed version {version}\n||",
             version=screed.__version__)

    log_info("|| If you use this script in a publication, please "
             "cite EACH of the following:\n||")

    if algorithm_list is None:
        algorithm_list = []

    algorithm_list.insert(0, 'software')

    for alg in algorithm_list:
        algstr = "||   * " + _algorithms[alg].encode(
            'utf-8', 'surrogateescape').decode('utf-8', 'replace')
        try:
            log_info(algstr)
        except UnicodeEncodeError:
            log_info(algstr.encode(sys.getfilesystemencoding(), 'replace'))

    log_info("||\n|| Please see http://khmer.readthedocs.org/en/"
             "latest/citations.html for details.\n")

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
