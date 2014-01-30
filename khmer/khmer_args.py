#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import os
import argparse

DEFAULT_K = 32
DEFAULT_N_HT = 4
DEFAULT_MIN_HASHSIZE = 1e6

def build_hash_args(descr=None):

    parser = argparse.ArgumentParser(
        description=descr,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
    env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', DEFAULT_MIN_HASHSIZE)

    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
   
    parser.add_argument('--ksize', '-k', type=int, dest='ksize',
                        default=env_ksize,
                        help='k-mer size to use')
    parser.add_argument('--n_hashes', '-N', type=int, dest='n_hashes',
                        default=env_n_hashes,
                        help='number of hash tables to use')
    parser.add_argument('--hashsize', '-x', type=float, dest='min_hashsize',
                        default=env_hashsize,
                        help='lower bound on hashsize to use')

    return parser


def build_counting_args(descr=None):

    if descr is None:
        descr = 'Build & load a counting Bloom filter.'

    parser = build_hash_args(descr=descr)

    return parser


def build_hashbits_args(descr=None):

    if descr is None:
        descr = 'Build & load a Bloom filter.'

    parser = build_hash_args(descr=descr)

    return parser

# deprecated, should use add_loadhash_args for input table
def build_counting_multifile_args():
    parser = argparse.ArgumentParser(
        description='Use a counting Bloom filter.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_table')
    parser.add_argument('input_filenames', nargs='+')

    return parser

# add an argument for loadhash with warning about parameters
def add_loadhash_args(parser):

    class LoadAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
            env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
            env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', 
                                           DEFAULT_MIN_HASHSIZE)

            from khmer.utils import print_error

            if getattr(namespace, 'ksize') != env_ksize or \
            getattr(namespace, 'n_hashes') != env_n_hashes or \
            getattr(namespace, 'min_hashsize') != env_hashsize:
                if values:
                    print_error('''
** WARNING: You are loading a saved hashtable from
{hash}, but have set hashtable parameters. 
Your values for ksize, n_hashes, and hashsize \
will be ignored.'''.format(hash=values))

            setattr(namespace, self.dest, values)

 

    parser.add_argument('-l', '--loadhash', dest='loadhash', default=None,
                        help='load a precomputed hashtable from disk',
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
    print_error(" - n hashes =     {0} \t\t(-N)".format(args.n_hashes))
    print_error(
        " - min hashsize = {0:5.2g} \t(-x)".format(args.min_hashsize)
    )
    print_error("")
    if hashtype == 'counting':
        print_error(
            "Estimated memory usage is {0:.2g} bytes "
            "(n_hashes x min_hashsize)".format(args.n_hashes * args.min_hashsize)
        )
    elif hashtype == 'hashbits':
        print_error(
            "Estimated memory usage is {0:.2g} bytes "
            "(n_hashes x min_hashsize / 8)".format(args.n_hashes * args.min_hashsize / 8)
        )
    
    print_error("-" * 8)
 
    if DEFAULT_MIN_HASHSIZE == args.min_hashsize and \
    not hasattr(args, 'loadhash'):
        print_error(
            "** WARNING: hashsize is default!  "
            "You absodefly want to increase this!\n** "
            "Please read the docs!\n"
        )


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
