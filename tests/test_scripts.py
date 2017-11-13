# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014-2015, Michigan State University.
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
# pylint: disable=C0111,C0103,E1103,unused-variable,protected-access


import csv
import json
import sys
import os
import stat
import threading
import io

import pytest
from . import khmer_tst_utils as utils
import khmer
from khmer import Countgraph, SmallCountgraph, Nodegraph
import khmer.kfile
import screed


def teardown():
    utils.cleanup()


def test_check_space():
    # @CTB this probably belongs in a new test file, along with other
    # tests of the file.py module.
    khmer.kfile.check_space(
        ['', utils.get_test_data('test-abund-read-2.fa')], False)


def test_load_into_counting():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 94' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_smallcount():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '--small-count']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 83' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_quiet():
    script = 'load-into-counting.py'
    args = ['-q', '-x', '1e3', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert len(out) == 0
    assert len(err) == 0
    assert os.path.exists(outfile)


def test_load_into_counting_autoargs_0():
    script = 'load-into-counting.py'

    outfile = utils.get_temp_filename('table')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args = ['-U', '1e7', '--fp-rate', '0.08', outfile, infile]
    (status, out, err) = utils.runscript(script, args)

    assert os.path.exists(outfile)
    assert 'INFO: Overriding default fp 0.1 with new fp: 0.08' in err, err
    assert ' tablesize is too small!' in err, err
    assert 'Estimated FP rate with current config is: 0.9999546' in err, err
    assert 'Recommended tablesize is: 1.77407e+07 bytes' in err, err


def test_load_into_counting_autoargs_1():
    script = 'load-into-counting.py'

    outfile = utils.get_temp_filename('table')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args = ['-U', '1e7', '--max-tablesize', '3e7', outfile, infile]
    (status, out, err) = utils.runscript(script, args)

    assert os.path.exists(outfile)
    assert "Ceiling is: 4.80833e+07 bytes" in err, err
    assert "set memory ceiling automatically." in err, err


def test_load_into_count_graphsize_warning():
    script = 'load-into-counting.py'
    args = ['-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert os.path.exists(outfile)
    assert "WARNING: tablesize is default!" in err


def test_load_into_counting_max_memory_usage_parameter():
    script = 'load-into-counting.py'
    args = ['-M', '2e3', '-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert os.path.exists(outfile)
    assert "WARNING: tablesize is default!" not in err

    kh = Countgraph.load(outfile)
    assert sum(kh.hashsizes()) < 3e8


def test_load_into_counting_abundance_dist_nobig():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20', '-b']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 94' in err, err
    assert os.path.exists(outfile)

    htfile = outfile
    outfile = utils.get_temp_filename('out')
    script2 = 'abundance-dist.py'
    args = ['-z', htfile, infile, outfile]
    (status, out, err) = utils.runscript(script2, args)
    assert 'WARNING: The loaded graph has bigcount' in err, err
    assert 'bigcount' in err, err


def test_load_into_counting_abundance_dist_squashing():
    graphfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args = [graphfile, infile]
    script = 'load-into-counting.py'
    utils.runscript(script, args)

    histogram = utils.get_temp_filename('histogram')
    infile = utils.get_test_data('test-abund-read-2.fa')
    args = [graphfile, infile, histogram]

    script = 'abundance-dist.py'
    # make histogram
    (status, out, err) = utils.runscript(script, args)
    assert os.path.exists(histogram)
    # attempt to overwrite histogram; fail
    failed = True
    try:
        (status, out, err) = utils.runscript(script, args)
        failed = False
    except AssertionError as error:
        assert "exists; not squashing" in str(error), str(error)

    assert failed, "Expected to fail"
    # attempt to overwrite with squashing; should work
    args = ['-s', graphfile, infile, histogram]
    (status, out, err) = utils.runscript(script, args)
    assert "squashing existing file" in err, err

    histfile = open(histogram, 'r')
    lines = histfile.readlines()
    # stripping because boo whitespace
    assert lines[1].strip() == "0,0,0,0.0", lines[1]
    assert lines[2].strip() == "1,83,83,1.0", lines[2]


# note: if run as root, will fail b/c root can write to anything
@pytest.mark.noroot
def test_load_into_counting_nonwritable():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('test-nonwritable')
    with open(outfile, 'w') as fout:
        fout.write("This file is non-writable (after this)")

    os.chmod(outfile, stat.S_IWOTH | stat.S_IRUSR)
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert 'does not have write permission; exiting' in err, err
    assert status == 1, status


@pytest.mark.huge
def test_load_into_counting_toobig():
    script = 'load-into-counting.py'
    args = ['-x', '1e12', '-N', '2', '-k', '20', '--force']

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == -1, status
    assert "MemoryError" in err, err


def test_load_into_counting_fail():
    script = 'load-into-counting.py'
    args = ['-x', '1e2', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    print(err)
    assert "** ERROR: the graph structure is too small" in err


def test_load_into_counting_multifile():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile, infile, infile, infile, infile,
                 infile, infile, infile, infile, infile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_tsv():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-s', 'tsv']

    outfile = utils.get_temp_filename('out.ct')
    tabfile = outfile + '.info.tsv'
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)
    assert os.path.exists(tabfile)
    with open(tabfile) as tabfh:
        tabfile_lines = tabfh.readlines()
    assert len(tabfile_lines) == 2
    outbase = os.path.basename(outfile)
    tsv = [outbase, '0.000', '95', '1001', infile]
    expected_tsv_line = '\t'.join(tsv) + '\n'
    assert tabfile_lines[1] == expected_tsv_line, tabfile_lines


def test_load_into_counting_json():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-s', 'json']

    outfile = utils.get_temp_filename('out.ct')
    jsonfile = outfile + '.info.json'
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)
    assert os.path.exists(jsonfile)

    with open(jsonfile) as jsonfh:
        got_json = json.load(jsonfh)
    outbase = os.path.basename(outfile)

    expected_json = {
        u"files": [infile],
        u"ht_name": outbase,
        u"num_kmers": 95,
        u"num_reads": 1001,
        u"fpr": 9.025048735197377e-11,
        u"mrinfo_version": "0.2.0",
    }

    assert got_json == expected_json, got_json


def test_load_into_counting_bad_summary_fmt():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-s', 'badfmt']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status != 0, status
    assert "invalid choice: 'badfmt'" in err, err


def test_load_into_counting_info_version():
    script = 'load-into-counting.py'
    args = ['-x', '1e5', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    info_file = outfile + '.info'
    assert os.path.exists(info_file), info_file
    with open(info_file) as info_fp:
        versionline = info_fp.readline()
    version = versionline.split(':')[1].strip()
    assert versionline.startswith('khmer version:'), versionline
    assert version == khmer.__version__, version


def _make_counting(infilename, SIZE=1e7, N=2, K=20, BIGCOUNT=True):
    script = 'load-into-counting.py'
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]

    if not BIGCOUNT:
        args.append('-b')

    outfile = utils.get_temp_filename('out.ct')

    args.extend([outfile, infilename])

    utils.runscript(script, args)
    assert os.path.exists(outfile)

    return outfile


def test_filter_stoptags():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)
    stopfile = utils.get_temp_filename('stoptags', in_dir)

    # first, copy test-abund-read-2.fa to 'test.fa' in the temp dir.
    # now, create a file with some stop tags in it --
    K = 18
    kh = khmer.Nodegraph(K, 1, 1)
    kh.add_stop_tag('GTTGACGGGGCTCAGGGG')
    kh.save_stop_tags(stopfile)
    del kh

    # finally, run filter-stoptags.
    script = 'filter-stoptags.py'
    args = ['-k', str(K), stopfile, infile, infile]
    utils.runscript(script, args, in_dir)

    # verify that the basic output file exists
    outfile = infile + '.stopfilt'
    assert os.path.exists(outfile), outfile

    # it should contain only one unique sequence, because we've trimmed
    # off everything after the beginning of the only long sequence in there.
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs


def test_filter_stoptags_fq():
    infile = utils.copy_test_data('test-abund-read-2.fq')
    in_dir = os.path.dirname(infile)
    stopfile = utils.get_temp_filename('stoptags', in_dir)

    # first, copy test-abund-read-2.fa to 'test.fa' in the temp dir.

    # now, create a file with some stop tags in it --
    K = 18
    kh = khmer.Nodegraph(K, 1, 1)
    kh.add_stop_tag('GTTGACGGGGCTCAGGGG')
    kh.save_stop_tags(stopfile)
    del kh

    # finally, run filter-stoptags.
    script = 'filter-stoptags.py'
    args = ['-k', str(K), stopfile, infile, infile]
    utils.runscript(script, args, in_dir)

    # verify that the basic output file exists
    outfile = infile + '.stopfilt'
    assert os.path.exists(outfile), outfile

    # it should contain only one unique sequence, because we've trimmed
    # off everything after the beginning of the only long sequence in there.
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs

    # make sure that record names are carried through unparsed
    names = [r.name for r in screed.open(outfile)]
    names = set(names)
    assert 'seq 1::BAR' in names


def test_count_median():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = infile + '.counts'

    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile).readlines()[1:]]
    data = set(data)
    assert len(data) == 2, data
    assert 'seq,1001,1001.0,0.0,18' in data, data
    assert '895:1:37:17593:9954/1,1,103.803741455,303.702941895,114' in data


def test_count_median_fq_csv():
    infile = utils.copy_test_data('test-abund-read-2.fq')
    outfile = infile + '.counts'

    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert len(data) == 4, data
    assert 'name,median,average,stddev,seqlen' in data
    assert 'seq,1001,1001.0,0.0,18' in data

    # verify that sequence names remain unparsed
    names = set([line.split(',')[0] for line in data])
    assert '895:1:37:17593:9954 1::FOO' in names, names


def test_count_median_fq_csv_stdout():
    infile = utils.copy_test_data('test-abund-read-2.fq')
    outfile = '-'

    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    (status, out, err) = utils.runscript(script, args)

    assert 'name,median,average,stddev,seqlen' in out
    assert 'seq,1001,1001.0,0.0,18' in out


def test_load_graph():
    script = 'load-graph.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    try:
        ht = Nodegraph.load(ht_file)
    except OSError as err:
        assert 0, str(err)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph.py worked :)
    subset = ht.do_subset_partition(0, 0)
    x = subset.count_partitions()
    assert x == (1, 0), x


@pytest.mark.known_failing
def test_oxli_build_graph():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = Nodegraph.load(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph.py worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


@pytest.mark.known_failing
def test_oxli_build_graph_unique_kmers_arg():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20', '-U', '3960']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err
    assert 'INFO: set memory ceiling automatically' in err, err
    assert 'Ceiling is: 1e+06 bytes' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = Nodegraph.load(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph.py worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


@pytest.mark.known_failing
def test_oxli_nocommand():
    script = 'oxli'

    (status, out, err) = utils.runscript(script, [])
    assert status == 0


def test_load_graph_no_tags():
    script = 'load-graph.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert Nodegraph.load(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...


@pytest.mark.known_failing
def test_oxli_build_graph_no_tags():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert Nodegraph.load(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...


def test_load_graph_fail():
    script = 'load-graph.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    assert "** ERROR: the graph structure is too small" in err


@pytest.mark.known_failing
def test_oxli_build_graph_fail():
    script = 'oxli'
    args = ['build-graph', '-x', '1e3', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    assert "** ERROR: the graph structure is too small" in err


@pytest.mark.known_failing
def test_oxli_build_graph_yuge():
    script = 'oxli'
    args = ['build-graph', '-M', '800T', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status != 0, status
    assert 'ERROR: Not enough free space on disk' in err


def test_load_graph_write_fp():
    script = 'load-graph.py'
    args = ['-x', '1e5', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    info_file = outfile + '.info'
    assert os.path.exists(info_file), info_file
    data = [x.strip() for x in open(info_file)]
    data = set(data)
    assert '3959 unique k-mers' in data, data
    assert 'false positive rate estimated to be 0.002' in data


@pytest.mark.known_failing
def test_oxli_build_graph_write_fp():
    script = 'oxli'
    # use small HT
    args = ['build-graph', '-x', '1e5', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    info_file = outfile + '.info'
    assert os.path.exists(info_file), info_file
    data = [x.strip() for x in open(info_file)]
    data = set(data)
    assert '3959 unique k-mers' in data
    assert 'false positive rate estimated to be 0.002' in data


def test_load_graph_multithread():
    script = 'load-graph.py'

    outfile = utils.get_temp_filename('test')
    infile = utils.get_test_data('test-reads.fa')

    args = ['-N', '4', '-x', '1e7', '-T', '8', outfile, infile]

    (status, out, err) = utils.runscript(script, args)


@pytest.mark.known_failing
def test_oxli_build_graph_multithread():
    script = 'oxli'

    outfile = utils.get_temp_filename('test')
    infile = utils.get_test_data('test-reads.fa')

    args = ['build-graph', '-N', '4', '-x', '1e7', '-T', '8', outfile, infile]

    (status, out, err) = utils.runscript(script, args)


def test_load_graph_max_memory_usage_parameter():
    script = 'load-graph.py'
    args = ['-M', '2e7', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    try:
        ht = Nodegraph.load(ht_file)
    except OSError as err:
        assert 0, str(err)

    assert (sum(ht.hashsizes()) / 8.) < 2e7, ht.hashsizes()


def _make_graph(infilename, min_hashsize=1e7, n_hashes=2, ksize=20,
                do_partition=False,
                annotate_partitions=False,
                stop_big_traverse=False):
    script = 'load-graph.py'
    args = ['-x', str(min_hashsize), '-N', str(n_hashes), '-k', str(ksize)]

    outfile = utils.get_temp_filename('out')
    infile = infilename

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    if do_partition:
        script = 'partition-graph.py'
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        utils.runscript(script, args)

        script = 'merge-partitions.py'
        args = [outfile, '-k', str(ksize)]
        utils.runscript(script, args)

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
            script = 'annotate-partitions.py'
            args = ["-k", str(ksize), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            utils.runscript(script, args, in_dir)

            baseinfile = os.path.basename(infilename)
            assert os.path.exists(os.path.join(in_dir, baseinfile + '.part'))

    return outfile


def test_partition_graph_1():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    utils.runscript('partition-graph.py', [graphbase])
    utils.runscript('merge-partitions.py', [graphbase, '-k', '20'])

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = Nodegraph.load(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_partition_graph_nojoin_k21():
    # test with K=21
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'), ksize=21)

    script = 'partition-graph.py'
    args = [graphbase]

    utils.runscript(script, args)

    script = 'merge-partitions.py'
    args = [graphbase, '-k', str(21)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = Nodegraph.load(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (99, 0), x          # should be 99 partitions at K=21


def test_partition_load_empty_pmap():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'), ksize=24)

    script = 'partition-graph.py'
    args = [graphbase, '-s', '10']

    utils.runscript(script, args)

    script = 'merge-partitions.py'
    args = [graphbase, '-k', '24']
    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == -1
    assert 'only a header and no partition IDs' in err


def test_partition_graph_nojoin_stoptags():
    # test with stoptags
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    # add in some stop tags
    ht = Nodegraph.load(graphbase)
    ht.add_stop_tag('TTGCATACGTTGAGCCAGCG')
    stoptags_file = graphbase + '.stoptags'
    ht.save_stop_tags(stoptags_file)
    del ht

    # run script with stoptags option
    script = 'partition-graph.py'
    args = ['--stoptags', stoptags_file, graphbase]

    utils.runscript(script, args)

    script = 'merge-partitions.py'
    args = [graphbase, '-k', str(20)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = Nodegraph.load(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (2, 0), x          # should be 2 partitions


def test_partition_graph_big_traverse():
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=False)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = Nodegraph.load(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_partition_graph_no_big_traverse():
    # do NOT exhaustively traverse
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=True)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = Nodegraph.load(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x[0] == 4, x       # should be four partitions, broken at knot.


def test_partition_find_knots_execute():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    script = 'partition-graph.py'
    args = [graphbase]

    utils.runscript(script, args)

    script = 'find-knots.py'
    args = [graphbase]
    utils.runscript(script, args)

    stoptags_file = graphbase + '.stoptags'
    assert os.path.exists(stoptags_file)


def test_partition_find_knots_existing_stoptags():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    script = 'partition-graph.py'
    args = [graphbase]
    utils.runscript(script, args)

    script = 'make-initial-stoptags.py'
    args = [graphbase]
    utils.runscript(script, args)

    script = 'find-knots.py'
    args = [graphbase]
    (status, out, err) = utils.runscript(script, args)

    stoptags_file = graphbase + '.stoptags'
    assert os.path.exists(stoptags_file)
    assert "loading stoptags" in err, err
    assert "these output stoptags will include the already" in err, err


def test_partition_graph_too_many_threads():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    utils.runscript('partition-graph.py', [graphbase, '--threads', '100'])
    utils.runscript('merge-partitions.py', [graphbase, '-k', '20'])

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = Nodegraph.load(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_annotate_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = 'annotate-partitions.py'
    args = ["-k", "20", graphbase, seqfile]
    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_annotate_partitions_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True,
                            ksize=21)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = 'annotate-partitions.py'
    args = ["-k", "21", graphbase, seqfile]
    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    print(parts)
    assert len(parts) == 99, len(parts)


def test_extract_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)


def test_extract_paired_inconsistent_formats():
    fa_seqfile = utils.get_test_data('random-20-a.fa')
    fq_seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        fa_seqfile, do_partition=True, annotate_partitions=True)
    fa_in_dir = os.path.dirname(graphbase)
    graphbase = _make_graph(
        fq_seqfile, do_partition=True, annotate_partitions=True)
    fq_in_dir = os.path.dirname(graphbase)
    # XXX
    # get the final part file
    fa_partfile = os.path.join(fa_in_dir, 'random-20-a.fa.part')
    fq_partfile = os.path.join(fq_in_dir, 'random-20-a.fq.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', fa_partfile, fq_partfile]

    failed = True
    try:
        utils.runscript(script, args, fa_in_dir)
        failed = False
    except AssertionError as err:
        assert "Input files must have consistent format." in str(err), err

    assert failed, "Expected to fail"


def test_extract_partitions_header_whitespace():
    seqfile = utils.get_test_data('test-overlap2.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'test-overlap2.fa.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '1 11960 11960 11960', dist.strip()

    parts = [r.name.split('\t')[1]
             for r in screed.open(partfile)]
    assert len(parts) == 13538, len(parts)
    parts = set(parts)
    assert len(parts) == 12602, len(parts)


def test_extract_partitions_fq():
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fq')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    screed_iter = screed.open(partfile)
    names = [r.name.split('\t')[0] for r in screed_iter]
    assert '35 1::FOO' in names
    assert '46 1::FIZ' in names

    screed_iter = screed.open(partfile)
    parts = [r.name.split('\t')[1] for r in screed_iter]

    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)

    quals = set([r.quality for r in screed.open(partfile)])
    quals = list(quals)
    assert quals[0], quals


def test_extract_partitions_output_unassigned():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-U', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    unassigned_file = os.path.join(in_dir, 'extracted.unassigned.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)
    assert os.path.exists(unassigned_file)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)


def test_extract_partitions_no_output_groups():
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-n', 'extracted', partfile]

    # We expect a sys.exit -> we need the test to be tolerant
    status, out, err = utils.runscript(script, args, in_dir)
    assert "NOT outputting groups! Beware!" in err
    # Group files are created after output_groups is
    # checked. They should not exist in this scenario
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert not os.path.exists(groupfile)


def test_extract_partitions_pid_0():
    partfile = utils.copy_test_data('random-20-a.fa.part')

    in_dir = os.path.dirname(partfile)
    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-U', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    unassigned_file = os.path.join(in_dir, 'extracted.unassigned.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)
    assert os.path.exists(unassigned_file)

    # Assert unassigned file not empty
    unassigned_content = open(unassigned_file).readline()
    assert unassigned_content.strip().split('\t')[0] != ''


def test_extract_partitions_multi_groups():
    partfile = utils.copy_test_data('random-20-a.fa.part')

    in_dir = os.path.dirname(partfile)

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-m', '1', '-X', '1', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    # Multiple group files are created after should be created
    groupfile1 = os.path.join(in_dir, 'extracted.group0000.fa')
    groupfile2 = os.path.join(in_dir, 'extracted.group0001.fa')
    groupfile3 = os.path.join(in_dir, 'extracted.group0002.fa')
    assert os.path.exists(groupfile1)
    assert os.path.exists(groupfile2)
    assert os.path.exists(groupfile3)


def test_extract_partitions_no_groups():
    empty_file = utils.copy_test_data('empty-file')
    in_dir = os.path.dirname(empty_file)

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', empty_file]

    status, _, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert "ERROR: Input file", "is empty; Exiting." in err
    assert status != 0
    # No group files should be created
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')

    assert not os.path.exists(groupfile)


def test_abundance_dist():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert (line == 'abundance,count,cumulative,cumulative_fraction'), line
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_quiet():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-z', '-q', htfile, infile, outfile]
    status, out, err = utils.runscript(script, args, in_dir)

    assert len(err) == 0

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert (line == 'abundance,count,cumulative,cumulative_fraction'), line
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_stdout():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-z', htfile, infile, "-"]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert '1,96,96,0.98' in out, out
    assert '1001,2,98,1.0' in out, out


def test_abundance_dist_nobigcount():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-b', '-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '255,2,98,1.0', line


def test_abundance_dist_threaded():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '--threads', '18',
            infile, outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_csv():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', infile,
            outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert (line == 'abundance,count,cumulative,cumulative_fraction'), line
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_nobigcount():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-b', infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '255,2,98,1.0', line


def test_abundance_dist_single_smallcount():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '--small-count',
            infile, outfile]
    utils.runscript(script, args, in_dir)


def test_abundance_dist_single_nosquash():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_quiet():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-q', '-x', '1e7', '-N', '2', '-k', '17', '-z', infile, outfile]
    status, out, err = utils.runscript(script, args, in_dir)

    assert len(err) == 0

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_savegraph():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    outfile = utils.get_temp_filename('test.dist')
    tabfile = utils.get_temp_filename('test-savegraph.ct')
    in_dir = os.path.dirname(infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '--savegraph',
            tabfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_do_partition():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "20", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_do_partition_no_big_traverse():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "20", "--no-big-traverse", "--threads=100", graphbase,
            seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_do_partition_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "21", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)

    assert len(parts) == 99, len(parts)


def test_do_partition_2_fq():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "21", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    screed_iter = screed.open(partfile)
    names = [r.name.split('\t')[0] for r in screed_iter]
    assert '35 1::FOO' in names
    assert '46 1::FIZ' in names


def test_make_initial_stoptags():
    # gen input files using load-graph.py -t
    # should keep test_data directory size down
    # or something like that
    # this assumes (obv.) load-graph.py works properly
    bzinfile = utils.copy_test_data('test-reads.fq.bz2')
    in_dir = os.path.dirname(bzinfile)

    genscript = 'load-graph.py'
    genscriptargs = ['test-reads', 'test-reads.fq.bz2']
    utils.runscript(genscript, genscriptargs, in_dir)

    # test input file gen'd by load-graph.pys
    infile = utils.get_temp_filename('test-reads.pt')
    infile2 = utils.get_temp_filename('test-reads.tagset', in_dir)

    # get file to compare against
    ex_outfile = utils.get_test_data('test-reads.stoptags')

    # actual output file
    outfile1 = utils.get_temp_filename('test-reads.stoptags', in_dir)

    script = 'make-initial-stoptags.py'
    # make-initial-stoptags has weird file argument syntax
    # read the code before modifying
    args = ['test-reads']

    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile1), outfile1


def test_make_initial_stoptags_load_stoptags():
    # gen input files using load-graph.py -t
    # should keep test_data directory size down
    # or something like that
    # this assumes (obv.) load-graph.py works properly
    bzinfile = utils.copy_test_data('test-reads.fq.bz2')
    in_dir = os.path.dirname(bzinfile)

    genscript = 'load-graph.py'
    genscriptargs = ['test-reads', 'test-reads.fq.bz2']
    utils.runscript(genscript, genscriptargs, in_dir)

    # test input file gen'd by load-graph.pys
    infile = utils.get_temp_filename('test-reads.pt')
    infile2 = utils.get_temp_filename('test-reads.tagset', in_dir)

    # get file to compare against
    ex_outfile = utils.get_test_data('test-reads.stoptags')

    # actual output file
    outfile1 = utils.get_temp_filename('test-reads.stoptags', in_dir)

    script = 'make-initial-stoptags.py'
    # make-initial-stoptags has weird file argument syntax
    # read the code before modifying
    args = ['test-reads']

    utils.runscript(script, args, in_dir)
    args = ['test-reads', '--stoptags', 'test-reads.stoptags']
    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile1), outfile1


def test_sample_reads_randomly():
    infile = utils.copy_test_data('test-reads.fa')
    in_dir = os.path.dirname(infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'850:2:1:1859:11742/1', '850:2:1:1859:11742/2',
                  '850:2:1:2131:17360/1', '850:2:1:2131:17360/2',
                  '850:2:1:2416:7565/1', '850:2:1:2416:7565/2',
                  '850:2:1:2490:13491/1', '850:2:1:2490:13491/2',
                  '850:2:1:2962:3999/1', '850:2:1:2962:3999/2',
                  '850:2:1:3096:20321/1', '850:2:1:3096:20321/2',
                  '850:2:1:3164:6414/1', '850:2:1:3164:6414/2',
                  '850:2:1:3206:13876/1', '850:2:1:3206:13876/2',
                  '850:2:1:3631:20919/1', '850:2:1:3631:20919/2',
                  '850:2:1:3655:15581/1', '850:2:1:3655:15581/2'}
    else:
        answer = {'850:2:1:1257:3404/1', '850:2:1:1257:3404/2',
                  '850:2:1:1362:19357/1', '850:2:1:1362:19357/2',
                  '850:2:1:1396:5659/1', '850:2:1:1396:5659/2',
                  '850:2:1:2063:11124/1', '850:2:1:2063:11124/2',
                  '850:2:1:2121:12070/1', '850:2:1:2121:12070/2',
                  '850:2:1:2528:15779/1', '850:2:1:2528:15779/2',
                  '850:2:1:2581:12886/1', '850:2:1:2581:12886/2',
                  '850:2:1:2864:8505/1', '850:2:1:2864:8505/2',
                  '850:2:1:3000:2015/1', '850:2:1:3000:2015/2',
                  '850:2:1:3302:5025/1', '850:2:1:3302:5025/2'}

    assert seqs == answer


def test_sample_reads_randomly_force_single():
    infile = utils.copy_test_data('test-reads.fa')
    in_dir = os.path.dirname(infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1', '--force_single']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'850:2:1:2399:20086/2',
                  '850:2:1:2273:13309/1',
                  '850:2:1:2065:16816/1',
                  '850:2:1:1984:7162/2',
                  '850:2:1:2691:14602/1',
                  '850:2:1:1762:5439/1',
                  '850:2:1:2503:4494/2',
                  '850:2:1:2263:11143/2',
                  '850:2:1:1792:15774/2',
                  '850:2:1:2084:17145/1'}
    else:
        answer = {'850:2:1:1199:4197/1',
                  '850:2:1:1251:16575/2',
                  '850:2:1:1267:6790/2',
                  '850:2:1:1601:4443/1',
                  '850:2:1:1625:19325/1',
                  '850:2:1:1832:14607/2',
                  '850:2:1:1946:20852/2',
                  '850:2:1:2401:4896/2',
                  '850:2:1:2562:1308/1',
                  '850:2:1:3123:15968/2'}

    assert seqs == answer


def test_sample_reads_randomly_force_single_outfile():
    infile = utils.copy_test_data('test-reads.fa')
    in_dir = os.path.dirname(infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1', '--force_single', '-o',
            in_dir + '/randreads.out']

    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = in_dir + '/randreads.out'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'850:2:1:2399:20086/2',
                  '850:2:1:2273:13309/1',
                  '850:2:1:2065:16816/1',
                  '850:2:1:1984:7162/2',
                  '850:2:1:2691:14602/1',
                  '850:2:1:1762:5439/1',
                  '850:2:1:2503:4494/2',
                  '850:2:1:2263:11143/2',
                  '850:2:1:1792:15774/2',
                  '850:2:1:2084:17145/1'}
    else:
        answer = {'850:2:1:1199:4197/1',
                  '850:2:1:1251:16575/2',
                  '850:2:1:1267:6790/2',
                  '850:2:1:1601:4443/1',
                  '850:2:1:1625:19325/1',
                  '850:2:1:1832:14607/2',
                  '850:2:1:1946:20852/2',
                  '850:2:1:2401:4896/2',
                  '850:2:1:2562:1308/1',
                  '850:2:1:3123:15968/2'}

    assert seqs == answer


def test_sample_reads_randomly_fq():
    infile = utils.copy_test_data('test-reads.fq.gz')
    in_dir = os.path.dirname(infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    if sys.version_info.major == 2:
        answer = {'850:2:1:2399:20086/2',
                  '850:2:1:1762:5439 1::FOO',
                  '850:2:1:2065:16816/1',
                  '850:2:1:2263:11143/2',
                  '850:2:1:1792:15774/2',
                  '850:2:1:2691:14602/1',
                  '850:2:1:2503:4494 1::FOO',
                  '850:2:1:2084:17145/1',
                  '850:2:1:1984:7162 1::FOO',
                  '850:2:1:2273:13309 1::FOO'}
    else:
        answer = {'850:2:1:1199:4197 1::FOO',
                  '850:2:1:1251:16575/2',
                  '850:2:1:1267:6790/2',
                  '850:2:1:1601:4443 1::FOO',
                  '850:2:1:1625:1932 1::FOO1',
                  '850:2:1:1832:14607 1::FOO',
                  '850:2:1:1946:20852 1::FOO',
                  '850:2:1:2401:4896/2',
                  '850:2:1:2562:1308/1',
                  '850:2:1:3123:15968/2'}

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))
    assert seqs == answer


def test_sample_reads_randomly_stdin_no_out():
    script = 'sample-reads-randomly.py'
    args = ['-']

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename" in err, err


def test_fastq_to_fasta():

    script = 'fastq-to-fasta.py'
    clean_infile = utils.copy_test_data('test-fastq-reads.fq')
    n_infile = utils.copy_test_data('test-fastq-n-reads.fq')

    clean_outfile = clean_infile + '.keep.fa'
    n_outfile = n_infile + '.keep.fa'

    in_dir = os.path.dirname(clean_infile)
    in_dir_n = os.path.dirname(n_infile)

    args = [clean_infile, '-n', '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0, len(out.splitlines())
    assert "No lines dropped" in err

    names = [r.name for r in screed.open(clean_outfile)]
    assert '895:1:1:1246:14654 1:N:0:NNNNN' in names, names

    args = [n_infile, '-n', '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 0
    assert "No lines dropped" in err

    args = [clean_infile, '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0
    assert "0 lines dropped" in err

    args = [n_infile, '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 0, out
    assert "4 lines dropped" in err, err

    args = [clean_infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) > 0
    assert "0 lines dropped" in err

    args = [n_infile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) > 0
    assert "4 lines dropped" in err

    args = [clean_infile, '-o', clean_outfile, '--gzip']
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0
    assert "0 lines dropped" in err

    args = [clean_infile, '-o', clean_outfile, '--bzip']
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0
    assert "0 lines dropped" in err


def test_fastq_to_fasta_streaming_compressed_gzip():

    script = 'fastq-to-fasta.py'
    infile = utils.copy_test_data('test-reads.fq.gz')
    in_dir = os.path.dirname(infile)
    fifo = utils.get_temp_filename('fifo')
    copyfilepath = utils.get_temp_filename('copied.fa.gz', in_dir)

    # make a fifo to simulate streaming
    os.mkfifo(fifo)
    args = ['--gzip', '-o', fifo, infile]
    # FIFOs MUST BE OPENED FOR READING BEFORE THEY ARE WRITTEN TO
    # If this isn't done, they will BLOCK and things will hang.
    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    copyfile = io.open(copyfilepath, 'wb')
    fifofile = io.open(fifo, 'rb')

    # read binary to handle compressed files
    chunk = fifofile.read(8192)
    while len(chunk) > 0:
        copyfile.write(chunk)
        chunk = fifofile.read(8192)

    fifofile.close()
    thread.join()
    copyfile.close()

    # verify that the seqs are there and not broken
    f = screed.open(copyfilepath)
    count = 0
    for _ in f:
        count += 1

    assert count == 25000, count
    f.close()

    # verify we're looking at a gzipped file
    gzfile = io.open(file=copyfilepath, mode='rb', buffering=8192)
    magic = b"\x1f\x8b\x08"  # gzip magic signature
    file_start = gzfile.peek(len(magic))
    assert file_start[:3] == magic, file_start[:3]


def test_fastq_to_fasta_streaming_compressed_bzip():

    script = 'fastq-to-fasta.py'
    infile = utils.copy_test_data('test-reads.fq.gz')
    in_dir = os.path.dirname(infile)
    fifo = utils.get_temp_filename('fifo')
    copyfilepath = utils.get_temp_filename('copied.fa.bz', in_dir)

    # make a fifo to simulate streaming
    os.mkfifo(fifo)
    args = ['--bzip', '-o', fifo, infile]
    # FIFOs MUST BE OPENED FOR READING BEFORE THEY ARE WRITTEN TO
    # If this isn't done, they will BLOCK and things will hang.
    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    copyfile = io.open(copyfilepath, 'wb')
    fifofile = io.open(fifo, 'rb')

    # read binary to handle compressed files
    chunk = fifofile.read(8192)
    while len(chunk) > 0:
        copyfile.write(chunk)
        chunk = fifofile.read(8192)

    fifofile.close()
    thread.join()
    copyfile.close()

    # verify that the seqs are there and not broken
    f = screed.open(copyfilepath)
    count = 0
    for _ in f:
        count += 1

    assert count == 25000, count
    f.close()

    # verify we're looking at a gzipped file
    bzfile = io.open(file=copyfilepath, mode='rb', buffering=8192)
    magic = b"\x42\x5a\x68"  # bzip magic signature
    file_start = bzfile.peek(len(magic))
    assert file_start[:3] == magic, file_start[:3]


def test_extract_long_sequences_fa():

    script = 'extract-long-sequences.py'
    fa_infile = utils.copy_test_data('paired-mixed.fa')

    fa_outfile = fa_infile + '.keep.fa'

    in_dir_fa = os.path.dirname(fa_infile)

    args = [fa_infile, '-l', '10', '-o', fa_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_fa)

    countlines = sum(1 for line in open(fa_outfile))
    assert countlines == 22, countlines

    names = [r.name for r in screed.open(fa_outfile)]
    assert "895:1:37:17593:9954/1" in names
    assert "895:1:37:17593:9954/2" in names


def test_extract_long_sequences_fq():

    script = 'extract-long-sequences.py'
    fq_infile = utils.copy_test_data('paired-mixed.fq')

    fq_outfile = fq_infile + '.keep.fq'

    in_dir_fq = os.path.dirname(fq_infile)

    args = [fq_infile, '-l', '10', '-o', fq_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_fq)

    countlines = sum(1 for line in open(fq_outfile))
    assert countlines == 44, countlines

    names = [r.name for r in screed.open(fq_outfile)]
    assert "895:1:37:17593:9954 1::foo" in names
    assert "895:1:37:17593:9954 2::foo" in names


def test_sample_reads_randomly_S():
    infile = utils.copy_test_data('test-fastq-reads.fq')
    in_dir = os.path.dirname(infile)

    script = 'sample-reads-randomly.py'

    # fix random number seed for reproducibility
    args = ['-N', '10', '-R', '1', '-S', '3']

    badargs = list(args)
    badargs.extend(['-o', 'test', infile, infile])
    (status, out, err) = utils.runscript(script, badargs, in_dir, fail_ok=True)
    assert status == 1, (status, out, err)
    assert "Error: cannot specify -o with more than one sample" in err

    args.append(infile)

    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset.0'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile, parse_description=True)])
    print(list(sorted(seqs)))

    print(seqs)
    if sys.version_info.major == 2:
        answer = {'895:1:1:1303:14389', '895:1:1:1347:3237',
                  '895:1:1:1295:6189', '895:1:1:1308:20421',
                  '895:1:1:1320:11648', '895:1:1:1352:5369',
                  '895:1:1:1318:10532', '895:1:1:1363:11839',
                  '895:1:1:1355:13535', '895:1:1:1349:15165'}
    else:
        answer = {'895:1:1:1290:11501', '895:1:1:1303:14389',
                  '895:1:1:1307:4308', '895:1:1:1308:2539',
                  '895:1:1:1331:1766', '895:1:1:1333:2512',
                  '895:1:1:1347:3237', '895:1:1:1363:11839',
                  '895:1:1:1378:18986', '895:1:1:1383:3089'}

    assert seqs == answer

    outfile = infile + '.subset.1'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile, parse_description=True)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = set(['895:1:1:1303:14389', '895:1:1:1373:4848',
                      '895:1:1:1357:19736', '895:1:1:1347:3237',
                      '895:1:1:1338:7557', '895:1:1:1388:11093',
                      '895:1:1:1296:1784', '895:1:1:1290:11501',
                      '895:1:1:1355:13535', '895:1:1:1303:6251'])
    else:
        answer = {'895:1:1:1255:18861', '895:1:1:1276:16426',
                  '895:1:1:1303:6251', '895:1:1:1308:20421',
                  '895:1:1:1314:10430', '895:1:1:1351:14718',
                  '895:1:1:1355:13535', '895:1:1:1358:4953',
                  '895:1:1:1362:3983', '895:1:1:1363:9988'}
    assert seqs == answer

    seqs = set([r.name for r in screed.open(outfile, parse_description=True)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'895:1:1:1303:14389', '895:1:1:1373:4848',
                  '895:1:1:1357:19736', '895:1:1:1347:3237',
                  '895:1:1:1338:7557', '895:1:1:1388:11093',
                  '895:1:1:1296:1784', '895:1:1:1290:11501',
                  '895:1:1:1355:13535', '895:1:1:1303:6251'}

    else:
        answer = {'895:1:1:1362:3983', '895:1:1:1363:9988',
                  '895:1:1:1314:10430', '895:1:1:1255:18861',
                  '895:1:1:1308:20421', '895:1:1:1358:4953',
                  '895:1:1:1351:14718', '895:1:1:1303:6251',
                  '895:1:1:1276:16426', '895:1:1:1355:13535'}

    assert seqs == answer


def execute_streaming_diginorm(ifilename):
    '''Helper function for the matrix of streaming tests for read_parser
    using diginorm, i.e. uncompressed fasta, gzip fasta, bz2 fasta,
    uncompressed fastq, etc.
    This is not directly executed but is run by the tests themselves
    '''
    # Get temp filenames, etc.
    fifo = utils.get_temp_filename('fifo')
    in_dir = os.path.dirname(fifo)
    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '-o', 'outfile', fifo]

    # make a fifo to simulate streaming
    os.mkfifo(fifo)

    # FIFOs MUST BE OPENED FOR READING BEFORE THEY ARE WRITTEN TO
    # If this isn't done, they will BLOCK and things will hang.
    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    ifile = io.open(ifilename, 'rb')
    fifofile = io.open(fifo, 'wb')
    # read binary to handle compressed files
    chunk = ifile.read(8192)
    while len(chunk) > 0:
        fifofile.write(chunk)
        chunk = ifile.read(8192)

    fifofile.close()

    thread.join()

    return in_dir + '/outfile'


def _execute_load_graph_streaming(filename):
    '''Helper function for the matrix of streaming tests using screed via
    filter-abund-single, i.e. uncompressed fasta, gzip fasta, bz2 fasta,
    uncompressed fastq, etc.
    This is not directly executed but is run by the tests themselves
    '''

    scripts = utils.scriptpath()
    infile = utils.copy_test_data(filename)
    in_dir = os.path.dirname(infile)

    args = '-x 1e7 -N 2 -k 20 out -'

    cmd = 'cat {infile} | {scripts}/load-graph.py {args}'.format(
        infile=infile, scripts=scripts, args=args)

    (status, out, err) = utils.run_shell_cmd(cmd, in_directory=in_dir)

    if status != 0:
        print(out)
        print(err)
        assert status == 0, status

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = os.path.join(in_dir, 'out')
    assert os.path.exists(ht_file), ht_file

    tagset_file = os.path.join(in_dir, 'out.tagset')
    assert os.path.exists(tagset_file), tagset_file

    ht = Nodegraph.load(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph.py worked :)
    subset = ht.do_subset_partition(0, 0)
    x = subset.count_partitions()
    assert x == (1, 0), x


def test_screed_streaming_ufa():
    # uncompressed fa
    o = execute_streaming_diginorm(utils.get_test_data('test-abund-read-2.fa'))

    pathstat = os.stat(o)
    seqs = [r.sequence for r in screed.open(o)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG')


def test_screed_streaming_ufq():
    # uncompressed fq
    o = execute_streaming_diginorm(utils.get_test_data('test-fastq-reads.fq'))

    seqs = [r.sequence for r in screed.open(o)]
    assert seqs[0].startswith('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')


def test_screed_streaming_bzipfq():
    # bzip compressed fq
    o = execute_streaming_diginorm(utils.get_test_data('100-reads.fq.bz2'))
    seqs = [r.sequence for r in screed.open(o)]
    assert len(seqs) == 100, seqs
    assert seqs[0].startswith('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'), seqs


def test_screed_streaming_bzipfa():
    # bzip compressed fa
    o = execute_streaming_diginorm(
        utils.get_test_data('test-abund-read-2.fa.bz2'))

    seqs = [r.sequence for r in screed.open(o)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG')


@pytest.mark.known_failing
def test_screed_streaming_gzipfq():
    # gzip compressed fq
    o = execute_streaming_diginorm(utils.get_test_data('100-reads.fq.gz'))
    assert os.path.exists(o)
    seqs = [r.sequence for r in screed.open(o)]
    assert seqs[0].startswith('CAGGCGCCCACCACCGTGCCCTCCAACCTG')


@pytest.mark.known_failing
def test_screed_streaming_gzipfa():
    o = execute_streaming_diginorm(
        utils.get_test_data('test-abund-read-2.fa.gz'))
    assert os.path.exists(o)
    seqs = [r.sequence for r in screed.open(o)]
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGG')


def test_read_parser_streaming_ufa():
    # uncompressed FASTA
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fa'))


def test_read_parser_streaming_ufq():
    # uncompressed FASTQ
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fq'))


@pytest.mark.known_failing
def test_read_parser_streaming_bzfq():
    # bzip compressed FASTQ
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fq.bz2'))


def test_read_parser_streaming_gzfq():
    # gzip compressed FASTQ
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fq.gz'))


@pytest.mark.known_failing
def test_read_parser_streaming_bzfa():
    # bzip compressed FASTA
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fa.bz2'))


def test_read_parser_streaming_gzfa():
    # gzip compressed FASTA
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fa.gz'))


def test_readstats():
    readstats_output = ("358 bp / 5 seqs; 71.6 average length",
                        "916 bp / 11 seqs; 83.3 average length")

    args = [utils.get_test_data("test-sweep-reads.fq"),
            utils.get_test_data("paired-mixed.fq")]
    status, out, err = utils.runscript('readstats.py', args)
    assert status == 0

    for k in readstats_output:
        assert k in out, (k, out)


def test_readstats_csv():
    readstats_output = ("358,5,71.6," +
                        utils.get_test_data("test-sweep-reads.fq"),
                        "916,11,83.3," +
                        utils.get_test_data("paired-mixed.fq"))

    args = [utils.get_test_data("test-sweep-reads.fq"),
            utils.get_test_data("paired-mixed.fq"),
            '--csv']
    status, out, err = utils.runscript('readstats.py', args)
    assert status == 0

    for k in readstats_output:
        assert k in out, (k, out)


def test_readstats_output():
    readstats_output = ("358 bp / 5 seqs; 71.6 average length",
                        "916 bp / 11 seqs; 83.3 average length")

    outfile = utils.get_temp_filename('output.txt')
    args = ["-o", outfile,
            utils.get_test_data("test-sweep-reads.fq"),
            utils.get_test_data("paired-mixed.fq")]

    status, _, _ = utils.runscript('readstats.py', args)
    assert status == 0

    out = open(outfile).read()

    for k in readstats_output:
        assert k in out, (k, out)


def test_readstats_empty():
    expected_output = "No sequences found in 2 files"

    args = [utils.get_test_data("test-empty.fa"),
            utils.get_test_data("test-empty.fa.bz2")]

    status, out, err = utils.runscript('readstats.py', args)
    assert status == 0

    assert expected_output in out


def test_trim_low_abund_1():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_trim_low_abund_1_duplicate_filename_err():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1', infile, infile]
    (status, out, err) = utils.runscript('trim-low-abund.py', args, in_dir,
                                         fail_ok=True)
    assert status == 1
    assert "Error: Cannot input the same filename multiple times." in str(err)


def test_trim_low_abund_1_stdin_err():
    args = ["-"]

    (status, out, err) = utils.runscript('trim-low-abund.py', args,
                                         fail_ok=True)
    assert status == 1
    assert "Accepting input from stdin; output filename must be provided" \
           in str(err)


def test_trim_low_abund_2():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    infile2 = utils.copy_test_data('test-abund-read-2.fa', 'copyDataTwo')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1', infile, infile2]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_trim_low_abund_2_o_gzip():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    infile2 = utils.copy_test_data('test-abund-read-2.fa', 'copyDataTwo')
    outfile = utils.get_temp_filename('out.gz')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1',
            "-o", outfile, "--gzip",
            infile, infile2]
    utils.runscript('trim-low-abund.py', args, in_dir)

    assert os.path.exists(outfile), outfile
    x = list(screed.open(outfile))
    assert len(x)

# make sure that FASTQ records are retained.


def test_trim_low_abund_3_fq_retained():
    infile = utils.copy_test_data('test-abund-read-2.fq')
    infile2 = utils.copy_test_data('test-abund-read-2.fq', 'copyDataTwo')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1', infile, infile2]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

    # check for 'quality' string.
    seqs = set([r.quality for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert '##################' in seqs


# test that the -V option does not trim sequences that are low abundance


def test_trim_low_abund_4_retain_low_abund():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# test that the -V option *does* trim sequences that are low abundance


def test_trim_low_abund_5_trim_high_abund():
    infile = utils.copy_test_data('test-abund-read-3.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # trimmed sequence @ error
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGC' in seqs

# test that -V/-Z setting - should not trip if -Z is set high enough.


def test_trim_low_abund_6_trim_high_abund_Z():
    infile = utils.copy_test_data('test-abund-read-3.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', '-Z', '25', infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # untrimmed seq.
    badseq = 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCgtgCCGCAGCTGTCGTCAGGG' \
             'GATTTCCGGGCGG'
    assert badseq in seqs       # should be there, untrimmed


def test_trim_low_abund_keep_paired():
    infile = utils.copy_test_data('test-abund-read-2.paired.fq')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-V", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = [r.name for r in screed.open(outfile)]
    assert seqs[-2:] == ['pair/1', 'pair/2'], seqs


def test_trim_low_abund_keep_paired_casava18():
    infile = utils.copy_test_data('test-abund-read-2.paired2.fq')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-V", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = [r.name for r in screed.open(outfile)]
    assert seqs[-2:] == ['pair:foo 1::N', 'pair:foo 2::N'], seqs


def test_trim_low_abund_highfpr():
    infile = utils.copy_test_data('test-abund-read-2.paired.fq')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1", "-N", "1", "-V", infile]
    code, out, err = utils.runscript('trim-low-abund.py', args, in_dir,
                                     fail_ok=True)

    assert code == 1
    assert '** ERROR: the graph structure is too small' in err, err


def test_trim_low_abund_trimtest():
    infile = utils.copy_test_data('test-abund-read-2.paired.fq')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-Z", "2", "-C", "1",
            "-V", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    for record in screed.open(outfile):
        if record.name == 'seqtrim/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCC'
        elif record.name == 'seqtrim/2':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGC'
        elif record.name == 'seqtrim2/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCA'


def test_trim_low_abund_trimtest_after_load():
    infile = utils.copy_test_data('test-abund-read-2.paired.fq')
    in_dir = os.path.dirname(infile)

    saved_table = utils.get_temp_filename('save.ct')

    args = ["-k", "17", "-x", "1e7", "-N", "2", saved_table, infile]
    utils.runscript('load-into-counting.py', args, in_dir)

    args = ["-Z", "2", "-C", "2", "-V", '--loadgraph', saved_table, infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    for record in screed.open(outfile):
        if record.name == 'seqtrim/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCC'
        elif record.name == 'seqtrim/2':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGC'
        elif record.name == 'seqtrim2/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCA'


def test_trim_low_abund_trimtest_savegraph():
    infile = utils.copy_test_data('test-abund-read-2.paired.fq')
    in_dir = os.path.dirname(infile)

    saved_table = utils.get_temp_filename('save.ct')

    args = ["-k", "17", "-x", "1e7", "-N", "2",
            "-Z", "2", "-C", "2", "-V", '--savegraph', saved_table, infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile
    assert os.path.exists(saved_table)

    for record in screed.open(outfile):
        if record.name == 'seqtrim/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCC'
        elif record.name == 'seqtrim/2':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGC'
        elif record.name == 'seqtrim2/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCA'


def test_trim_low_abund_no_summary_info_by_default():
    infile = utils.copy_test_data("test-abund-read-2.fa")
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-o", "summary", infile]
    _, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    summary_fname = os.path.join(in_dir, "summary.info.json")
    print(os.path.exists(summary_fname))
    assert not os.path.exists(summary_fname), summary_fname


def test_trim_low_abund_summary_info_json():
    # test JSON file with summary info is created
    infile = utils.copy_test_data("test-abund-read-2.fa")
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "--summary-info", "json",
            "-o", "summary", infile]
    _, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    summary_fname = os.path.join(in_dir, "summary.info.json")
    assert os.path.exists(summary_fname), summary_fname
    with open(summary_fname) as f:
        assert json.load(f), 'summary file does not contain valid JSON'


def test_trim_low_abund_summary_info_tsv():
    # test TSV file with summary info is created
    infile = utils.copy_test_data("test-abund-read-2.fa")
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "--summary-info", "tsv",
            "-o", "summary", infile]
    _, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    summary_fname = os.path.join(in_dir, "summary.info.tsv")
    assert os.path.exists(summary_fname), summary_fname
    with open(summary_fname) as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        lines = [row for row in reader]
        assert len(lines) == 1


# test that -o/--out option outputs to STDOUT
def test_trim_low_abund_stdout():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", infile, "-o", "-"]
    _, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    # attempt to parse output to check it is in FASTA format
    stream = io.StringIO(out)
    assert list(screed.fasta.fasta_iter(stream)), "can't parse stdout"

    # can't test that the correct message appears because we redirect
    # the output when under testing. Instead check that incorrect message
    # does not appear.
    assert 'output in *.abundtrim' not in err


def test_trim_low_abund_output_named():
    # check the output filename is mentioned when it is explicitly set
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", infile,
            "-o", "explicitname.abundtrim"]
    _, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    assert 'output in explicitname.abundtrim' in err


def test_trim_low_abund_diginorm_coverage_err():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-M", "1e7", infile, "--diginorm-coverage", "21"]
    status, out, err = utils.runscript('trim-low-abund.py', args, in_dir,
                                       fail_ok=True)

    print(out, err)
    assert status == 1
    assert 'Error: --diginorm-coverage given, but --diginorm not specified.' \
           in err, err


def test_trim_low_abund_diginorm_single_pass():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-M", "1e7", infile, "--diginorm", "--single-pass"]
    status, out, err = utils.runscript('trim-low-abund.py', args, in_dir,
                                       fail_ok=True)

    assert status == 1
    assert "Error: --diginorm and --single-pass are incompatible!" \
           in err, err


def test_trim_low_abund_varcov_err():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-M", "1e7", infile, "-Z", "21"]
    status, out, err = utils.runscript('trim-low-abund.py', args, in_dir,
                                       fail_ok=True)

    print(out, err)
    assert status == 1
    assert 'Error: --trim-at-coverage/-Z given' in err, err


def test_trim_low_abund_single_pass():
    infile = utils.copy_test_data('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    args = ["-M", "1e7", infile, "-V", '--single-pass']
    status, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    assert status == 0


def test_trim_low_abund_quiet():
    infile = utils.copy_test_data('test-reads.fa')
    in_dir = os.path.dirname(infile)

    args = ["-q", "-M", "1e7", infile, "-V", '-Z', '5', '-C', '1']
    status, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    assert status == 0
    assert len(out) == 0
    assert len(err) == 0


def test_trim_low_abund_reporting():
    infile = utils.copy_test_data('test-reads.fa')
    in_dir = os.path.dirname(infile)

    args = ["-M", "1e7", infile, "-V", '-Z', '5', '-C', '1']
    status, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    assert status == 0
    assert '11157 11161 848236 2 152' in err


def test_roundtrip_casava_format_1():
    # check to make sure that extract-paired-reads produces a file identical
    # to the input file when only paired data is given.

    infile = utils.copy_test_data('casava_18-pe.fq')
    in_dir = os.path.dirname(infile)

    _, out, err = utils.runscript('extract-paired-reads.py', [infile], in_dir)

    r = open(infile).read()

    outfile = infile + '.pe'
    r2 = open(outfile).read()
    assert r == r2, (r, r2)


def test_roundtrip_casava_format_2():
    # check that split-paired-reads -> interleave-reads produces a file
    # identical to input, when only paired reads are given.

    infile = utils.copy_test_data('casava_18-pe.fq')
    outfile = utils.get_temp_filename('test2.fq')
    in_dir = os.path.dirname(infile)

    _, out, err = utils.runscript('split-paired-reads.py', [infile], in_dir)

    utils.runscript('interleave-reads.py', [infile + '.1',
                                            infile + '.2',
                                            '-o', outfile], in_dir)

    r = open(infile).read()
    r2 = open(outfile).read()
    assert r == r2, (r, r2)


def test_existence_failure():
    expected_output = 'ERROR: Input file'

    args = [utils.get_temp_filename('thisfiledoesnotexistatall')]

    status, out, err = utils.runscript(
        'extract-paired-reads.py', args, fail_ok=True)
    assert status == 1

    assert expected_output in err


def test_roundtrip_commented_format():
    """Split/interleave roundtrip for old style format with comments (#873).

    This should produce a file identical to the input when only paired
    reads are given.
    """
    infile = utils.copy_test_data('old-style-format-w-comments.fq')
    outfile = utils.get_temp_filename('test2.fq')
    in_dir = os.path.dirname(infile)

    _, out, err = utils.runscript('split-paired-reads.py', [infile], in_dir)

    utils.runscript('interleave-reads.py', [infile + '.1',
                                            infile + '.2',
                                            '-o', outfile], in_dir)

    r = open(infile).read()
    r2 = open(outfile).read()
    assert r == r2, (r, r2)


def test_unique_kmers_defaults():
    infile = utils.copy_test_data('random-20-a.fa')

    args = ['-k', '20', '-e', '0.01', infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'.format(infile)
            in err)
    assert 'Total estimated number of unique 20-mers: 3950' in err


def test_unique_kmers_report_fp():
    infile = utils.copy_test_data('random-20-a.fa')
    outfile = utils.get_temp_filename('report.unique')

    args = ['-k', '20', '-e', '0.01', '-R', outfile, infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'.format(infile)
            in err)
    assert 'Total estimated number of unique 20-mers: 3950' in err

    with open(outfile, 'r') as report_fp:
        outf = report_fp.read().splitlines()
        assert '3950 20 (total)' in outf
        assert '3950 20 total' in outf


def test_unique_kmers_diagnostics():
    infile = utils.copy_test_data('random-20-a.fa')

    args = ['-k', '20', '-e', '0.01', '--diagnostics', infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    out = out.splitlines()
    assert ('expected_fp\tnumber_hashtable(Z)\t'
            'size_hashtable(H)\texpected_memory_usage' in err)


def test_unique_kmers_multiple_inputs():
    infiles = []
    for fname in ('random-20-a.fa', 'paired-mixed.fa'):
        infile = utils.copy_test_data(fname)
        infiles.append(infile)

    args = ['-k', '20', '-e', '0.01']
    args += infiles

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'
            .format(infiles[0]) in err)
    assert ('Estimated number of unique 20-mers in {0}: 232'.format(infiles[1])
            in err)
    assert 'Total estimated number of unique 20-mers: 4170' in err


@pytest.mark.parametrize("scriptname",
                         [entry for entry in os.listdir(utils.scriptpath())
                          if entry.endswith('.py')])
def test_version_and_basic_citation(scriptname):
    with open(os.path.join(utils.scriptpath(), scriptname)) as script:
        line = script.readline()
        line = script.readline()
        if 'khmer' in line:
            # check citation information appears when using --info
            status, out, err = utils.runscript(scriptname, ["--info"])
            assert status == 0, status
            print(out)
            print(err)
            assert "publication" in err, err
            assert "usage:" not in err, err

            # check citation information appears in --version
            status, out, err = utils.runscript(scriptname, ["--version"])
            assert status == 0, status
            assert "publication" in err, err
            assert "usage:" not in err, err

            # check citation information appears in --help
            status, out, err = utils.runscript(scriptname, ["--help"])
            assert status == 0, status
            assert "publication" in err, err
            assert "usage:" in out, out
