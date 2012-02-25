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

def test_load_into_counting_fail():
    script = scriptpath('load-into-counting.py')
    args = ['-x', '1e2', '-N', '2', '-k', '20'] # use small HT
    
    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = runscript(script, args)
    print out
    print err
    assert status == -1
    assert "ERROR:" in err

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

    outfile = infile + '.minkeep'
    assert os.path.exists(outfile), outfile

    seqs = [ r.sequence for r in screed.open(outfile) ]
    assert len(seqs) == 5, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs

def test_count_median():
    infile = utils.get_temp_filename('test.fa')
    outfile = infile + '.counts'
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = scriptpath('count-median.py')
    args = [counting_ht, infile, outfile]
    (status, out, err) = runscript(script, args)
    assert status == 0

    assert os.path.exists(outfile), outfile

    data = [ x.strip() for x in open(outfile) ]
    data = set(data)
    assert len(data) == 2, data
    assert 'seq 1001 1001.0 0.0 18' in data
    assert '895:1:37:17593:9954/1 1 103.803741455 303.702941895 114' in data

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

def test_load_graph_fail():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e3', '-N', '2', '-k', '20'] # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = runscript(script, args)
    assert status == -1
    assert "ERROR:" in err

def _make_graph(infilename, SIZE=1e7, N=2, K=20,
                do_partition=False,
                annotate_partitions=False,
                stop_big_traverse=False):
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

    if do_partition:
        script = scriptpath('partition-graph.py')
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        (status, out, err) = runscript(script, args)
        print out
        print err
        assert status == 0

        script = scriptpath('merge-partitions.py')
        args = [outfile, '-k', str(K)]
        (status, out, err) = runscript(script, args)
        print out
        print err
        assert status == 0

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
            script = scriptpath('annotate-partitions.py')
            args = ["-k", str(K), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            (status, out, err) = runscript(script, args, in_dir)
            assert status == 0

            baseinfile = os.path.basename(infilename)
            assert os.path.exists(os.path.join(in_dir, baseinfile + '.part'))

    return outfile

def test_partition_graph_1():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('partition-graph.py')
    args = [graphbase]

    (status, out, err) = runscript(script, args)
    assert status == 0

    script = scriptpath('merge-partitions.py')
    args = [graphbase, '-k', str(20)]
    (status, out, err) = runscript(script, args)
    print out
    print err
    assert status == 0

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0)          # should be exactly one partition.

def test_partition_graph_nojoin_k21():
    # test with K=21
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'), K=21)
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('partition-graph.py')
    args = [graphbase]

    (status, out, err) = runscript(script, args)
    assert status == 0

    script = scriptpath('merge-partitions.py')
    args = [graphbase, '-k', str(21)]
    (status, out, err) = runscript(script, args)
    print out
    print err
    assert status == 0

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (99, 0)          # should be 99 partitions at K=21

def test_partition_graph_nojoin_stoptags():
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

    script = scriptpath('merge-partitions.py')
    args = [graphbase, '-k', str(20)]
    (status, out, err) = runscript(script, args)
    print out
    print err
    assert status == 0

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (2, 0)          # should be 2 partitions

def test_partition_graph_big_traverse():
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=False)
    in_dir = os.path.dirname(graphbase)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0)          # should be exactly one partition.

def test_partition_graph_no_big_traverse():
    # do NOT exhaustively traverse
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=True)
    in_dir = os.path.dirname(graphbase)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.ht')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (4, 0), x       # should be four partitions, broken at knot.

def test_annotate_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = scriptpath('annotate-partitions.py')
    args = ["-k", "20", graphbase, seqfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [ r.name.split('\t')[1] for r in screed.open(partfile) ]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1

def test_annotate_partitions_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True,
                            K=21)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = scriptpath('annotate-partitions.py')
    args = ["-k", "21", graphbase, seqfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [ r.name.split('\t')[1] for r in screed.open(partfile) ]
    parts = set(parts)
    print parts
    assert len(parts) == 99, len(parts)

def test_extract_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['extracted', partfile]
    
    (status, out, err) = runscript(script, args, in_dir)
    print out
    print err
    assert status == 0

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [ r.name.split('\t')[1] for r in screed.open(partfile) ]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)

def test_abundance_dist():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17)

    script = scriptpath('abundance-dist.py')
    args = ['-z', htfile, infile, outfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '1001 2 98 1.0', line

def test_count_overlap_with_curve():
    seqfile1 = utils.get_test_data('test-overlap1.fa')
    seqfile2 = utils.get_test_data('test-overlap2.fa')
    in_dir = os.path.dirname(seqfile1)
    script = scriptpath('count-overlap.py')
    curvefile = seqfile1+'.curve'
    outfile = seqfile1+'.out'
    args = ['--ksize', '32', '--n_hashes', '4', '--hashsize','2000000000',\
            '--curve', curvefile,seqfile1,seqfile2,outfile]
    (status, out, err) = runscript(script, args, in_dir)
    print out, err
    assert status == 0
    assert os.path.exists(outfile), outfile
    assert os.path.exists(curvefile),curvefile
    #report file
    data = [ x.strip() for x in open(outfile) ]
    data = set(data)
    assert len(data) == 10, data
    assert '# of unique k-mers: 440346' in data
    assert '# of occupied bin: 440299' in data
    assert '# of unique k-mers: 581866' in data
    assert '# of occupied bin: 581783' in data
    assert 'false positive rate: 7.160103e-15' in data
    assert '# of overlap unique k-mers: 184849' in data
    #curve file
    data = [ x.strip() for x in open(curvefile) ]
    data = set(data)
    assert len(data) == 100, data
    assert '6021 0' in data
    assert '29649 40' in data
    assert '471277 74260' in data
    assert '529993 132976' in data
    assert '581866 184849' in data

def test_count_overlap_without_curve():
    seqfile1 = utils.get_test_data('test-overlap1.fa')
    seqfile2 = utils.get_test_data('test-overlap2.fa')
    in_dir = os.path.dirname(seqfile1)
    script = scriptpath('count-overlap.py')
    outfile = seqfile1+'.out'
    args = ['--ksize', '32', '--n_hashes', '4', '--hashsize','2000000000',\
            seqfile1,seqfile2,outfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0
    assert os.path.exists(outfile), outfile
    #report file
    data = [ x.strip() for x in open(outfile) ]
    data = set(data)
    assert len(data) == 10, data
    assert '# of unique k-mers: 440346' in data
    assert '# of occupied bin: 440299' in data
    assert '# of unique k-mers: 581866' in data
    assert '# of occupied bin: 581783' in data
    assert 'false positive rate: 7.160103e-15' in data
    assert '# of overlap unique k-mers: 184849' in data

def test_count_overlap_cpp():
    seqfile1 = utils.get_test_data('test-overlap1.ht')
    seqfile2 = utils.get_test_data('test-overlap2.fa')
    in_dir = os.path.dirname(seqfile1)
    script = scriptpath('count-overlap_cpp.py')
    outfile = seqfile1+'.out'
    args = ['--ksize', '32', '--n_hashes', '4', '--hashsize','2000000000',\
            seqfile1,seqfile2,outfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0
    assert os.path.exists(outfile), outfile
    #report file
    data = [ x.strip() for x in open(outfile) ]
    data = set(data)
    assert 'unique k-mers in dataset2:581866' in data
    assert 'overlap k-mers:184849' in data

