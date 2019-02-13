# vim: set encoding=utf-8
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
"""Common argparse constructs."""


import sys
import argparse
import math
import textwrap
from argparse import _VersionAction
from collections import namedtuple
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import screed
import khmer
from khmer import extract_countgraph_info
from khmer import __version__
from .utils import print_error
from .khmer_logger import log_info, log_warn, configure_logging


DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = 1e6
DEFAULT_N_THREADS = 1

ALGORITHMS = {
    'software': 'MR Crusoe et al., '
                '2015. https://doi.org/10.12688/f1000research.6924.1',
    'diginorm': 'CT Brown et al., arXiv:1203.4802 [q-bio.GN]',
    'streaming': 'Q Zhang, S Awad, CT Brown, '
                 'https://doi.org/10.7287/peerj.preprints.890v1',
    'graph': 'J Pell et al., https://doi.org/10.1073/pnas.1121464109',
    'counting': 'Q Zhang et al., '
                'https://doi.org/10.1371/journal.pone.0101271',
    'sweep': 'C Scott, MR Crusoe, and CT Brown, unpublished',
    'SeqAn': 'A. Döring et al. https://doi.org:80/10.1186/1471-2105-9-11',
    'hll': 'Irber and Brown. https://doi.org/10.1101/056846'
}


class CitationAction(argparse.Action):
    # pylint: disable=too-few-public-methods
    """Output citation information and exit."""

    def __init__(self, *args, **kwargs):
        self.citations = kwargs.pop('citations')
        super(CitationAction, self).__init__(*args, nargs=0,
                                             default=argparse.SUPPRESS,
                                             **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        info(parser.prog, self.citations)
        parser.exit()


class _HelpAction(argparse._HelpAction):
    # pylint: disable=too-few-public-methods, protected-access
    def __call__(self, parser, namespace, values, option_string=None):
        info(parser.prog, parser._citations)
        super(_HelpAction, self).__call__(parser, namespace, values,
                                          option_string=option_string)


class _VersionStdErrAction(_VersionAction):
    # pylint: disable=too-few-public-methods, protected-access
    """Force output to StdErr."""

    def __call__(self, parser, namespace, values, option_string=None):
        # have to call info() directly as the version action exits
        # which means parse_args() does not get a chance to run
        info(parser.prog, parser._citations)
        version = self.version
        if version is None:
            version = parser.version
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stderr)
        parser.exit()


class ComboFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    """Both ArgumentDefaults and RawDescription formatters."""

    pass


class KhmerArgumentParser(argparse.ArgumentParser):
    """Specialize ArgumentParser with khmer defaults.

    Take care of common arguments and setup printing of citation information.
    """

    def __init__(self, citations=None, formatter_class=ComboFormatter,
                 **kwargs):
        super(KhmerArgumentParser, self).__init__(
            formatter_class=formatter_class, add_help=False, **kwargs)
        self._citations = citations

        self.add_argument('--version', action=_VersionStdErrAction,
                          version='khmer {v}'.format(v=__version__))
        self.add_argument('--info', action=CitationAction,
                          help='print citation information',
                          citations=self._citations)
        self.add_argument('-h', '--help', action=_HelpAction,
                          default=argparse.SUPPRESS,
                          help='show this help message and exit')

    def parse_args(self, args=None, namespace=None):
        args = super(KhmerArgumentParser, self).parse_args(args=args,
                                                           namespace=namespace)

        # some scripts do not have a quiet flag, assume quiet=False for those
        if 'quiet' not in args or not args.quiet:
            info(self.prog, self._citations)

        return args


# Temporary fix to argparse FileType which ignores the
# binary mode flag. Upstream bug tracked in https://bugs.python.org/issue14156
# pylint: disable=too-few-public-methods,missing-docstring
class FileType(argparse.FileType):
    def __call__(self, fname):
        # detect if stdout is being faked (StringIO during unit tests) in
        # which case we do not have to do anything
        if (fname == '-' and
                sys.version_info.major == 3 and
                not isinstance(sys.stdout, StringIO)):
            if 'r' in self._mode:
                fname = sys.stdin.fileno()
            elif 'w' in self._mode:
                fname = sys.stdout.fileno()

        return super(FileType, self).__call__(fname)


def memory_setting(label):
    """
    Parse user-supplied memory setting.

    Converts strings into floats, representing a number of bytes. Supports the
    following notations.
      - raw integers: 1, 1000, 1000000000
      - scientific notation: 1, 1e3, 1e9
      - "common" notation: 1, 1K, 1G

    Suffixes supported: K/k, M/m, G/g, T/t. Do not include a trailing B/b.
    """
    suffixes = {
        'K': 1000.0,
        'M': 1000.0 ** 2,
        'G': 1000.0 ** 3,
        'T': 1000.0 ** 4,
    }
    try:
        mem = float(label)
        return mem
    except ValueError:
        prefix = label[:-1]
        suffix = label[-1:].upper()
        if suffix not in suffixes.keys():
            raise ValueError('cannot parse memory setting "{}"'.format(label))
        try:
            multiplier = float(prefix)
            return multiplier * suffixes[suffix]
        except ValueError:
            raise ValueError('cannot parse memory setting "{}"'.format(label))


def optimal_size(num_kmers, mem_cap=None, fp_rate=None):
    """
    Utility function for estimating optimal countgraph args.

      - num_kmers: number of unique kmers [required]
      - mem_cap: the allotted amount of memory [optional, conflicts with f]
      - fp_rate: the desired false positive rate [optional, conflicts with M]
    """
    if all((num_kmers is not None, mem_cap is not None, fp_rate is None)):
        return estimate_optimal_with_K_and_M(num_kmers, mem_cap)
    elif all((num_kmers is not None, mem_cap is None, fp_rate is not None)):
        return estimate_optimal_with_K_and_f(num_kmers, fp_rate)
    else:
        raise TypeError("num_kmers and either mem_cap or fp_rate"
                        " must be defined.")


def check_conflicting_args(args, hashtype):
    """
    Check argparse args object for conflicts.

    e.g. --loadgraph and --ksize being set.
    """
    if getattr(args, "quiet", None):
        configure_logging(args.quiet)

    loadgraph_table_conflicts = {"ksize": DEFAULT_K,
                                 "n_tables": DEFAULT_N_TABLES,
                                 "max_tablesize": DEFAULT_MAX_TABLESIZE}

    loadgraph_autoarg_conflicts = ("unique_kmers", "max_memory_usage")

    if getattr(args, "loadgraph", None):

        # check for table config args
        for key, value in loadgraph_table_conflicts.items():
            if getattr(args, key, value) != value:
                log_warn('''
*** WARNING: You are loading a saved k-mer countgraph from
*** {hashfile}, but have set k-mer table parameters.
*** Your values for ksize, n_tables, and tablesize
*** will be ignored.'''.format(hashfile=args.loadgraph))
                break  # no repeat warnings

        for element in loadgraph_autoarg_conflicts:
            if getattr(args, element, None):
                log_warn("\n*** WARNING: You have asked that the graph size be"
                         " automatically calculated\n"
                         "*** (by using -U or -M).\n"
                         "*** But you are loading an existing graph!\n"
                         "*** Size will NOT be set automatically.")
                break  # no repeat warnings

        infoset = None
        if hashtype in ('countgraph', 'smallcountgraph'):
            infoset = extract_countgraph_info(args.loadgraph)
        if infoset is not None:
            ksize = infoset.ksize
            max_tablesize = infoset.table_size
            n_tables = infoset.n_tables
            args.ksize = ksize
            args.n_tables = n_tables
            args.max_tablesize = max_tablesize
            if infoset.ht_type == khmer.FILETYPES['SMALLCOUNT']:
                args.small_count = True


def check_argument_range(low, high, parameter_name):
    """Check if parameter value is in the range `low` to `high`."""
    def _in_range(value):
        value = int(value)
        if not low <= value < high:
            print_error("\n** ERROR: khmer only supports "
                        "%i <= %s < %i.\n" % (low, parameter_name, high))
            sys.exit(1)

        else:
            return value

    return _in_range


# pylint: disable=invalid-name
def estimate_optimal_with_K_and_M(num_kmers, mem_cap):
    """
    Estimate optimal countgraph args.

     - num_kmers: number of unique kmer
     - mem_cap: the allotted amount of memory
    """
    n_tables = math.log(2) * (mem_cap / float(num_kmers))
    int_n_tables = int(n_tables)
    if int_n_tables == 0:
        int_n_tables = 1
    ht_size = int(mem_cap / int_n_tables)
    mem_cap = ht_size * int_n_tables
    fp_rate = (1 - math.exp(-num_kmers / float(ht_size))) ** int_n_tables
    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(int_n_tables, ht_size, mem_cap, fp_rate)


# pylint: disable=invalid-name
def estimate_optimal_with_K_and_f(num_kmers, des_fp_rate):
    """
    Estimate optimal memory.

    - num_kmers: the number of unique kmers
    - des_fp_rate: the desired false positive rate
    """
    n_tables = math.log(des_fp_rate, 0.5)
    int_n_tables = int(n_tables)
    if int_n_tables == 0:
        int_n_tables = 1

    ht_size = int(-num_kmers / (
        math.log(1 - des_fp_rate ** (1 / float(int_n_tables)))))
    mem_cap = ht_size * int_n_tables
    fp_rate = (1 - math.exp(-num_kmers / float(ht_size))) ** int_n_tables

    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(int_n_tables, ht_size, mem_cap, fp_rate)


def graphsize_args_report(unique_kmers, fp_rate):
    """
    Assemble output string for optimal arg sandbox scripts.

    - unique_kmers: number of uniqe k-mers
    - fp_rate: desired false positive rate
    """
    to_print = []

    to_print.append('')  # blank line
    to_print.append('number of unique k-mers: \t{0}'.format(unique_kmers))
    to_print.append('false positive rate: \t{:>.3f}'.format(fp_rate))
    to_print.append('')  # blank line
    to_print.append('If you have expected false positive rate to achieve:')
    to_print.append('expected_fp\tnumber_hashtable(Z)\tsize_hashtable(H)\t'
                    'expected_memory_usage')

    for fp_rate in range(1, 10):
        num_tables, table_size, mem_cap, fp_rate = \
            optimal_size(unique_kmers, fp_rate=fp_rate / 10.0)
        to_print.append('{:11.3f}\t{:19}\t{:17e}\t{:21e}'.format(fp_rate,
                                                                 num_tables,
                                                                 table_size,
                                                                 mem_cap))

    mem_list = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000, 5000]

    to_print.append('')  # blank line
    to_print.append('If you have expected memory to use:')
    to_print.append('expected_memory_usage\tnumber_hashtable(Z)\t'
                    'size_hashtable(H)\texpected_fp')

    for mem in mem_list:
        num_tables, table_size, mem_cap, fp_rate =\
            optimal_size(unique_kmers, mem_cap=mem * 1000000000)
        to_print.append('{:21e}\t{:19}\t{:17e}\t{:11.3f}'.format(mem_cap,
                                                                 num_tables,
                                                                 table_size,
                                                                 fp_rate))
    return "\n".join(to_print)


def _check_fp_rate(args, desired_max_fp):
    """
    Check if the desired_max_fp rate makes sense.

    - args: argparse args object, possible members: fp_rate, max_memory_usage,
      unique_kmers, force, max_tablesize
    - desired_max_fp: desired maximum false positive rate
    """
    if not args.unique_kmers:
        return args

    # Do overriding of default script FP rate
    if args.fp_rate:
        log_info("*** INFO: Overriding default fp {def_fp} with new fp:"
                 " {new_fp}", def_fp=desired_max_fp, new_fp=args.fp_rate)
        desired_max_fp = args.fp_rate

    # If we have the info we need to work with, do the stuff
    if args.max_memory_usage:
        # verify that this is a sane memory usage restriction
        res = estimate_optimal_with_K_and_M(args.unique_kmers,
                                            args.max_memory_usage)
        if res.fp_rate > desired_max_fp:
            print("""
*** ERROR: The given restrictions yield an estimate false positive rate of {0},
*** which is above the recommended false positive ceiling of {1}!"""
                  .format(res.fp_rate, desired_max_fp), file=sys.stderr)
            if not args.force:
                print("NOTE: This can be overridden using the --force"
                      " argument", file=sys.stderr)
                print("*** Aborting...!", file=sys.stderr)
                sys.exit(1)
    else:
        res = estimate_optimal_with_K_and_f(args.unique_kmers,
                                            desired_max_fp)
        if args.max_tablesize and args.max_tablesize < res.htable_size:
            log_warn("\n*** Warning: The given tablesize is too small!")
            log_warn("*** Recommended tablesize is: {tsize:5g} bytes",
                     tsize=res.htable_size)
            log_warn("*** Current is: {tsize:5g} bytes",
                     tsize=args.max_tablesize)
            res = estimate_optimal_with_K_and_M(args.unique_kmers,
                                                args.max_tablesize)
            log_warn("*** Estimated FP rate with current config is: {fp}\n",
                     fp=res.fp_rate)
        else:
            if res.mem_use < 1e6:  # one megabyteish
                args.max_memory_usage = 1e6
            else:
                args.max_memory_usage = res.mem_use
            log_info("*** INFO: set memory ceiling automatically.")
            log_info("*** Ceiling is: {ceil:3g} bytes\n",
                     ceil=float(args.max_memory_usage))
            args.max_mem = res.mem_use

    return args


def build_graph_args(descr=None, epilog=None, parser=None, citations=None):
    """Build an ArgumentParser with args for bloom filter based scripts."""
    expert_help = '--help-expert' in sys.argv
    if expert_help:
        sys.argv.append('--help')

    if parser is None:
        parser = KhmerArgumentParser(description=descr, epilog=epilog,
                                     citations=citations)

    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K,
                        help='k-mer size to use')

    help = ('number of tables to use in k-mer countgraph' if expert_help
            else argparse.SUPPRESS)
    parser.add_argument('--n_tables', '-N', type=int,
                        default=DEFAULT_N_TABLES,
                        help=help)

    parser.add_argument('-U', '--unique-kmers', type=float, default=0,
                        help='approximate number of unique kmers in the input'
                             ' set')
    parser.add_argument('--fp-rate', type=float, default=None,
                        help="Override the automatic FP rate setting for the"
                        " current script")

    group = parser.add_mutually_exclusive_group()
    help = ('upper bound on tablesize to use; overrides --max-memory-usage/-M'
            if expert_help else argparse.SUPPRESS)
    group.add_argument('--max-tablesize', '-x', type=float,
                       default=DEFAULT_MAX_TABLESIZE,
                       help=help)
    group.add_argument('-M', '--max-memory-usage', type=memory_setting,
                       help='maximum amount of memory to use for data ' +
                       'structure')

    return parser


def build_counting_args(descr=None, epilog=None, citations=None):
    """Build an ArgumentParser with args for countgraph based scripts."""
    parser = build_graph_args(descr=descr, epilog=epilog, citations=citations)

    parser.add_argument('--small-count', default=False, action='store_true',
                        help='Reduce memory usage by using a smaller counter'
                        ' for individual kmers.')

    return parser


def build_nodegraph_args(descr=None, epilog=None, parser=None, citations=None):
    """Build an ArgumentParser with args for nodegraph based scripts."""
    parser = build_graph_args(descr=descr, epilog=epilog, parser=parser,
                              citations=citations)

    return parser


def add_loadgraph_args(parser):
    """Common loadgraph argument."""
    parser.add_argument('-l', '--loadgraph', metavar="filename", default=None,
                        help='load a precomputed k-mer graph from disk')


def calculate_graphsize(args, graphtype, multiplier=1.0):
    """
    Transform the table parameters into a size.

    The return value refers to the target size (in buckets, not bytes) of each
    individual table in the graph.
    """
    if graphtype not in khmer._buckets_per_byte:
        raise ValueError('unknown graph type: ' + graphtype)

    if args.max_memory_usage:
        tablesize = float(multiplier) * (khmer._buckets_per_byte[graphtype] *
                                         args.max_memory_usage / args.n_tables)
    else:
        tablesize = args.max_tablesize

    return tablesize


def create_nodegraph(args, ksize=None, multiplier=1.0, fp_rate=0.01):
    """Create and return a nodegraph."""
    args = _check_fp_rate(args, fp_rate)

    if hasattr(args, 'force'):
        if args.n_tables > 20:
            if not args.force:
                print_error(
                    "\n** ERROR: khmer only supports number "
                    "of tables <= 20.\n")
                sys.exit(1)
            else:
                log_warn("\n*** Warning: Maximum recommended number of "
                         "tables is 20, discarded by force nonetheless!\n")

    if ksize is None:
        ksize = args.ksize
    if ksize > 32:
        print_error("\n** ERROR: khmer only supports k-mer sizes <= 32.\n")
        sys.exit(1)

    tablesize = calculate_graphsize(args, 'nodegraph', multiplier)
    return khmer.Nodegraph(ksize, tablesize, args.n_tables)


def create_countgraph(args, ksize=None, multiplier=1.0, fp_rate=0.1):
    """Create and return a countgraph."""
    args = _check_fp_rate(args, fp_rate)

    if hasattr(args, 'force'):
        if args.n_tables > 20:
            if not args.force:
                print_error(
                    "\n** ERROR: khmer only supports number "
                    "of tables <= 20.\n")
                sys.exit(1)
            else:
                if args.n_tables > 20:
                    log_warn("\n*** Warning: Maximum recommended number of "
                             "tables is 20, discarded by force nonetheless!\n")

    if ksize is None:
        ksize = args.ksize
    if ksize > 32:
        print_error("\n** ERROR: khmer only supports k-mer sizes <= 32.\n")
        sys.exit(1)

    if args.small_count:
        tablesize = calculate_graphsize(args, 'smallcountgraph',
                                        multiplier=multiplier)
        return khmer.SmallCountgraph(ksize, tablesize, args.n_tables)
    else:
        tablesize = calculate_graphsize(args, 'countgraph',
                                        multiplier=multiplier)
        cg = khmer.Countgraph(ksize, tablesize, args.n_tables)
        if hasattr(args, 'bigcount'):
            cg.set_use_bigcount(args.bigcount)
        return cg


def create_matching_nodegraph(countgraph):
    """Create a Nodegraph matched in size to a Countgraph

    Use this to create Nodegraphs for kmer tracking and similar. The
    created Nodegraph will have the same number of buckets in its
    tables as `countgraph`.
    """
    tablesizes = countgraph.hashsizes()
    return khmer.Nodegraph(countgraph.ksize(), 1, 1, primes=tablesizes)


def report_on_config(args, graphtype='countgraph'):
    """Print out configuration.

    Summarize the configuration produced by the command-line arguments
    made available by this module.
    """
    check_conflicting_args(args, graphtype)
    if graphtype not in khmer._buckets_per_byte:
        raise ValueError('unknown graph type: ' + graphtype)

    tablesize = calculate_graphsize(args, graphtype)
    maxmem = args.n_tables * tablesize / khmer._buckets_per_byte[graphtype]
    log_info("\nPARAMETERS:")
    log_info(" - kmer size =     {ksize} \t\t(-k)", ksize=args.ksize)
    log_info(" - n tables =      {ntables} \t\t(-N)", ntables=args.n_tables)
    log_info(" - max tablesize = {tsize:5.2g} \t(-x)", tsize=tablesize)
    log_info("Estimated memory usage is {mem:.1f} Gb "
             "({bytes:.2g} bytes = {ntables} bytes x {tsize:5.2g} entries "
             "/ {div:d} entries per byte)", bytes=maxmem, mem=maxmem / 1e9,
             div=khmer._buckets_per_byte[graphtype], ntables=args.n_tables,
             tsize=tablesize)
    log_info("-" * 8)

    if tablesize == DEFAULT_MAX_TABLESIZE and \
       not getattr(args, 'loadgraph', None):
        log_warn('''\

** WARNING: tablesize is default!
** You probably want to increase this with -M/--max-memory-usage!
** Please read the docs!
''')


def add_threading_args(parser):
    """Add option for threading to options parser."""
    parser.add_argument('-T', '--threads', default=DEFAULT_N_THREADS, type=int,
                        help='Number of simultaneous threads to execute')


def sanitize_help(parser):
    """Remove Sphinx directives & reflow text to width of 79 characters."""
    wrapper = textwrap.TextWrapper(width=79)
    parser.description = wrapper.fill(parser.description)
    if not parser.epilog:
        return parser
    cleanlog = parser.epilog.replace(':option:', '').replace(
        ':program:', '').replace('::', ':').replace('``', '"')
    newlog = prev_section = ""
    for section in cleanlog.split('\n\n'):
        if section.startswith('    '):
            newlog += section + '\n'
        else:
            if prev_section.startswith('    '):
                newlog += '\n'
            newlog += wrapper.fill(section) + '\n\n'
        prev_section = section
    parser.epilog = newlog
    return parser


def info(scriptname, algorithm_list=None):
    """Print version and project info to stderr."""
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
        algstr = "||   * " + ALGORITHMS[alg].encode(
            'utf-8', 'surrogateescape').decode('utf-8', 'replace')
        try:
            log_info(algstr)
        except UnicodeEncodeError:
            log_info(algstr.encode(sys.getfilesystemencoding(), 'replace'))

    log_info("||\n|| Please see http://khmer.readthedocs.io/en/"
             "latest/citations.html for details.\n")

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
