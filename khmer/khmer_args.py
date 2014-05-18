'''
This file is part of khmer, http://github.com/ged-lab/khmer/, and is
Copyright (C) Michigan State University, 2009-2014. It is licensed under
the three-clause BSD license; see doc/LICENSE.txt.
Contact: khmer-project@idyll.org
'''
import sys
import os
import argparse
from khmer import extract_countinghash_info, extract_hashbits_info
from khmer import __version__
import screed

DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MIN_TABLESIZE = 1e6
DEFAULT_N_THREADS = 1


class ComboFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    pass


def build_hash_args(descr=None, epilog=None):
    """Build an argparse.ArgumentParser with arguments for hash* based
    scripts and return it.
    """

    parser = argparse.ArgumentParser(
        description=descr, epilog=epilog,
        formatter_class=ComboFormatter)

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_tables = os.environ.get('KHMER_N_TABLES', DEFAULT_N_TABLES)
    env_tablesize = os.environ.get('KHMER_MIN_TABLESIZE',
                                   DEFAULT_MIN_TABLESIZE)

    parser.add_argument('--version', action='version',
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
    """Build an argparse.ArgumentParser with arguments for counting_hash
    based scripts and return it.
    """

    parser = build_hash_args(descr=descr, epilog=epilog)
    parser.hashtype = 'counting'

    return parser


def build_hashbits_args(descr=None, epilog=None):
    """Build an argparse.ArgumentParser with arguments for hashbits based
    scripts and return it.
    """

    parser = build_hash_args(descr=descr, epilog=epilog)
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
    """
        Summarizes the configuration produced by the command-line arguments
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
    parser.add_argument('--threads', '-T', default=DEFAULT_N_THREADS, type=int,
                        help='Number of simultaneous threads to execute')

_algorithms = {
    'software': 'MR Crusoe et al., 2014. doi: 10.6084/m9.figshare.979190',
    'diginorm': "CT Brown et al., arXiv:1203.4802 [q-bio.GN]",
    'graph': "J Pell et al., PNAS, 2014 (PMID 22847406)",
    'counting': "Q Zhang et al., arXiv:1309.2975 [q-bio.GN]",
    'sweep': 'C Scott, MR Crusoe, and CT Brown, unpublished'
}


def info(scriptname, algorithm_list=None):
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
        sys.stderr.write(_algorithms[alg])
        sys.stderr.write("\n")

    sys.stderr.write("||\n|| Please see the CITATION file for details.\n\n")

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
