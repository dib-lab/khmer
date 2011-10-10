import sys, os, shutil
from cStringIO import StringIO
import traceback

import khmer_tst_utils as utils
import khmer
import screed

def scriptpath(scriptfile):
    path = os.path.join(utils.thisdir, '..', 'scripts', scriptfile)
    path = os.path.abspath(path)
    return path

def teardown():
    utils.cleanup()

def runscript(scriptname, args, in_directory=None):
    """
    Run the given Python script, with the given args, in the given directory,
    using 'execfile'.
    """
    sysargs = [scriptname]
    sysargs.extend(args)

    cwd = os.getcwd()

    try:
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stderr = StringIO()

        if in_directory:
            os.chdir(in_directory)

        try:
            print 'running:', scriptname, 'in:', in_directory
            execfile(scriptname, { '__name__' : '__main__' })
            status = 0
        except:
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        os.chdir(cwd)

    return status, out, err

####

def test_load_into_counting():
    script = scriptpath('load-into-counting.py')
    args = ['-x', '1e7', '-N', '2', '-k', '20']
    
    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = runscript(script, args)
    assert status == 0
    assert os.path.exists(outfile)

def _make_counting(infilename, SIZE=1e7, N=2, K=20):
    script = scriptpath('load-into-counting.py')
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]
    
    outfile = utils.get_temp_filename('out.kh')

    args.extend([outfile, infilename])

    (status, out, err) = runscript(script, args)
    assert status == 0
    assert os.path.exists(outfile)

    return outfile

def test_filter_abund_1():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = [counting_ht, infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

def test_filter_abund_2():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-C', '1', counting_ht, infile, infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

def test_filter_stoptags():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    stopfile = utils.get_temp_filename('stoptags', in_dir)

    # first, copy test-abund-read-2.fa to 'test.fa' in the temp dir.
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    # now, create a file with some stop tags in it --
    K = 18
    kh = khmer.new_hashbits(K, 1, 1)
    kh.add_stop_tag('GTTGACGGGGCTCAGGGG')
    kh.save_stop_tags(stopfile)
    del kh
    
    # finally, run filter-stoptags.
    script = scriptpath('filter-stoptags.py')
    args = ['-k', str(K), stopfile, infile, infile]
    (status, out, err) = runscript(script, args, in_dir)
    print out
    print err
    assert status == 0

    # verify that the basic output file exists
    outfile = infile + '.stopfilt'
    assert os.path.exists(outfile), outfile

    # it should contain only one unique sequence, because we've trimmed
    # off everything after the beginning of the only long sequence in there.
    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs

def test_normalize_by_median():
    CUTOFF='1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [ r.sequence for r in screed.open(outfile) ]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs

def test_normalize_by_median_2():
    CUTOFF='2'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [ r.sequence for r in screed.open(outfile) ]
    assert len(seqs) == 2, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1] == 'GGTTGACGGGGCTCAGGG', seqs

def test_normalize_by_min():
    CUTOFF='5'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('normalize-by-min.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [ r.sequence for r in screed.open(outfile) ]
    assert len(seqs) == 5, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs

####

def test_load_graph():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = runscript(script, args)
    assert status == 0

    ht_file = outfile + '.ht'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_hashbits(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x

def _make_graph(infilename, SIZE=1e7, N=2, K=20):
    script = scriptpath('load-graph.py')
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data(infilename)

    args.extend([outfile, infile])

    (status, out, err) = runscript(script, args)
    assert status == 0

    ht_file = outfile + '.ht'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    return outfile

def test_partition_graph_1():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('partition-graph.py')
    args = [graphbase]

    (status, out, err) = runscript(script, args)
    assert status == 0

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0)          # should be exactly one partition.

def test_partition_graph_2():
    # test with K=21
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'), K=21)
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('partition-graph.py')
    args = [graphbase]

    (status, out, err) = runscript(script, args)
    assert status == 0

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (99, 0)          # should be 99 partitions at K=21

def test_partition_graph_3():
    # test with stoptags
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))
    in_dir = os.path.dirname(graphbase)

    # add in some stop tags
    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.add_stop_tag('TTGCATACGTTGAGCCAGCG')
    stoptags_file = graphbase + '.stoptags'
    ht.save_stop_tags(stoptags_file)
    del ht

    # run script with stoptags option
    script = scriptpath('partition-graph.py')
    args = ['--stoptags', stoptags_file, graphbase]

    (status, out, err) = runscript(script, args)
    assert status == 0

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (2, 0)          # should be 2 partitions
