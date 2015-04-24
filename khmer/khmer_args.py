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
from khmer import extract_countinghash_info, extract_hashbits_info
from khmer import __version__
import screed

DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MIN_TABLESIZE = 1e6
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


def build_hash_args(descr=None, epilog=None, parser=None):
    """Build an ArgumentParser with args for bloom filter based scripts."""
    if parser is None:
        parser = argparse.ArgumentParser(description=descr, epilog=epilog,
                                         formatter_class=ComboFormatter)

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_tables = os.environ.get('KHMER_N_TABLES', DEFAULT_N_TABLES)
    env_tablesize = os.environ.get('KHMER_MIN_TABLESIZE',
                                   DEFAULT_MIN_TABLESIZE)

    parser.add_argument('--version', action=_VersionStdErrAction,
                        version='khmer {v}'.format(v=__version__))
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')

    parser.add_argument('--ksize', '-k', type=int, default=env_ksize,
                        help='k-mer size to use')
    parser.add_argument('--n_tables', '-N', type=int,
                        default=env_n_tables,
                        help='number of k-mer counting tables to use')
    parser.add_argument('--min-tablesize', '-x', type=float,
                        default=env_tablesize,
                        help='lower bound on tablesize to use')

    return parser


def build_counting_args(descr=None, epilog=None):
    """Build an ArgumentParser with args for counting_hash based scripts."""
    parser = build_hash_args(descr=descr, epilog=epilog)
    parser.hashtype = 'counting'

    return parser


def build_hashbits_args(descr=None, epilog=None, parser=None):
    """Build an ArgumentParser with args for hashbits based scripts."""

    parser = build_hash_args(descr=descr, epilog=epilog, parser=parser)
    parser.hashtype = 'hashbits'

    return parser

# add an argument for loadhash with warning about parameters


def add_loadhash_args(parser):

    class LoadAction(argparse.Action):

        def __call__(self, parser, namespace, values, option_string=None):
            env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
            env_n_tables = os.environ.get('KHMER_N_TABLES', DEFAULT_N_TABLES)
            env_tablesize = os.environ.get('KHMER_MIN_TABLESIZE',
                                           DEFAULT_MIN_TABLESIZE)

            from khmer.utils import print_error

            setattr(namespace, self.dest, values)

            if getattr(namespace, 'ksize') != env_ksize or \
               getattr(namespace, 'n_tables') != env_n_tables or \
               getattr(namespace, 'min_tablesize') != env_tablesize:
                if values:
                    print_error('''
** WARNING: You are loading a saved k-mer table from
{hashfile}, but have set k-mer table parameters.
Your values for ksize, n_tables, and tablesize
will be ignored.'''.format(hashfile=values))

            if hasattr(parser, 'hashtype'):
                info = None
                if parser.hashtype == 'hashbits':
                    info = extract_hashbits_info(
                        getattr(namespace, self.dest))
                elif parser.hashtype == 'counting':
                    info = extract_countinghash_info(
                        getattr(namespace, self.dest))
                if info:
                    K = info[0]
                    x = info[1]
                    n = info[2]
                    setattr(namespace, 'ksize', K)
                    setattr(namespace, 'n_tables', n)
                    setattr(namespace, 'min_tablesize', x)

    parser.add_argument('-l', '--loadtable', metavar="filename", default=None,
                        help='load a precomputed k-mer table from disk',
                        action=LoadAction)


def report_on_config(args, hashtype='counting'):
    """Print out configuration.

    Summarize the configuration produced by the command-line arguments
    made available by this module.
    """
    from khmer.utils import print_error

    if args.quiet:
        return

    print_error("\nPARAMETERS:")
    print_error(" - kmer size =    {0} \t\t(-k)".format(args.ksize))
    print_error(" - n tables =     {0} \t\t(-N)".format(args.n_tables))
    print_error(
        " - min tablesize = {0:5.2g} \t(-x)".format(args.min_tablesize)
    )
    print_error("")
    if hashtype == 'counting':
        print_error(
            "Estimated memory usage is {0:.2g} bytes "
            "(n_tables x min_tablesize)".format(
                args.n_tables * args.min_tablesize))
    elif hashtype == 'hashbits':
        print_error(
            "Estimated memory usage is {0:.2g} bytes "
            "(n_tables x min_tablesize / 8)".format(args.n_tables *
                                                    args.min_tablesize / 8)
        )

    print_error("-" * 8)

    if DEFAULT_MIN_TABLESIZE == args.min_tablesize and \
       not hasattr(args, 'loadtable'):
        print_error(
            "** WARNING: tablesize is default!  "
            "You absodefly want to increase this!\n** "
            "Please read the docs!\n"
        )


def add_threading_args(parser):
    """Add option for threading to options parser."""
    parser.add_argument('--threads', '-T', default=DEFAULT_N_THREADS, type=int,
                        help='Number of simultaneous threads to execute')

_algorithms = {
    'software': 'MR Crusoe et al., '
    '2014. http://dx.doi.org/10.6084/m9.figshare.979190',
    'diginorm': 'CT Brown et al., arXiv:1203.4802 [q-bio.GN]',
    'streaming': 'Q Zhang, S Awad, CT Brown, unpublished',
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

    sys.stderr.write("\n")
    sys.stderr.write("|| This is the script '%s' in khmer.\n"
                     "|| You are running khmer version %s\n" %
                     (scriptname, khmer.__version__,))
    sys.stderr.write("|| You are also using screed version %s\n||\n"
                     % screed.__version__)

    sys.stderr.write("|| If you use this script in a publication, please "
                     "cite EACH of the following:\n||\n")

    if algorithm_list is None:
        algorithm_list = []

    algorithm_list.insert(0, 'software')

    for alg in algorithm_list:
        sys.stderr.write("||   * ")
        algstr = _algorithms[alg].encode(
            'utf-8', 'surrogateescape').decode('utf-8', 'replace')
        try:
            sys.stderr.write(algstr)
        except UnicodeEncodeError:
            sys.stderr.write(
                algstr.encode(sys.getfilesystemencoding(), 'replace'))
        sys.stderr.write("\n")

    sys.stderr.write("||\n|| Please see http://khmer.readthedocs.org/en/"
                     "latest/citations.html for details.\n\n")

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
