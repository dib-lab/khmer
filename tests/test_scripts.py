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
            print 'arguments', sysargs
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

def DEBUG_runscript(scriptname, args, in_directory=None):
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

        os.chdir(cwd)

    return status, "", ""

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

def _make_counting(infilename, SIZE=1e7, N=2, K=20, BIGCOUNT=True):
    script = scriptpath('load-into-counting.py')
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]

    if not BIGCOUNT:
        args.append('-b')
    
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

# make sure that FASTQ records are retained.
def test_filter_abund_3_fq_retained():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
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

    # check for 'accuracy' string.
    seqs = set([ r.accuracy for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs
    assert '##################' in seqs

def test_filter_abund_1_singlefile():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('filter-abund-single.py')
    args = ['-x', '1e7', '-N', '2', '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# test that the -V option does not trim sequences that are low abundance
def test_filter_abund_4_retain_low_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-V', counting_ht, infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# test that the -V option *does* trim sequences that are low abundance
def test_filter_abund_5_trim_high_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-V', counting_ht, infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs

    # trimmed sequence @ error
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGC' in seqs

# test that -V/-Z setting - should not trip if -Z is set high enough.
def test_filter_abund_6_trim_high_abund_Z():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-V', '-Z', '25', counting_ht, infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs

    # untrimmed seq.
    badseq = 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCgtgCCGCAGCTGTCGTCAGGGGATTTCCGGGCGG'
    assert badseq in seqs       # should be there, untrimmed

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

    seqs = sorted([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs
    
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs
    assert seqs[1].startswith('GGTTGACGGGGCTCAGGGGG'), seqs

def test_normalize_by_median_paired():
    CUTOFF='1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-paired.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    print out
    print err
    assert status == 0

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [ r.sequence for r in screed.open(outfile) ]
    assert len(seqs) == 2, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1].startswith('GGTTGACGGGGCTCAGGG'), seqs

def test_normalize_by_median_impaired():
    CUTOFF='1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    
    shutil.copyfile(utils.get_test_data('test-abund-read-impaired.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status != 0

def test_normalize_by_median_force():
    CUTOFF='1'
    
    corrupt_infile = utils.get_temp_filename('test-corrupt.fq')
    good_infile = utils.get_temp_filename('test-good.fq',
                                        tempdir=os.path.dirname(corrupt_infile))
    
    in_dir = os.path.dirname(good_infile)
    
    shutil.copyfile(utils.get_test_data('test-error-reads.fq'), corrupt_infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), good_infile)
    
    script = scriptpath('normalize-by-median.py')
    args = ['-f', '-C', CUTOFF, '-k', '17', corrupt_infile, good_infile]
    
    (status, out, err) = runscript(script, args, in_dir)

    assert os.path.exists(corrupt_infile + '.ht.failed')
    
    test_ht = khmer.load_counting_hash(corrupt_infile + '.ht.failed')
    test_good_read = 'CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'
    test_good_read2 = 'TAGTATCATCAAGGTTCAAGATGTTAATGAATAACAATTGCGCAGCAA'
    assert test_ht.count(test_good_read[:17]) > 0
    assert test_ht.count(test_good_read2[:17]) > 0
    assert status == 0
    assert os.path.exists(corrupt_infile + '.ht.failed')
    assert '*** Skipping' in err
    assert '** IOErrors' in err

def test_normalize_by_median_dumpfrequency():
    CUTOFF='1'
    
    infiles = [utils.get_temp_filename('test-0.fq')]
    in_dir = os.path.dirname(infiles[0])
    for x in range(1,5):
        infiles.append(utils.get_temp_filename('test-{}.fq'.format(x),
                                                tempdir=in_dir))
    
    for infile in infiles:
        shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)
    
    script = scriptpath('normalize-by-median.py')
    args = ['-d', '2', '-C', CUTOFF, '-k', '17']
    args.extend(infiles)
    
    (status, out, err) = runscript(script, args, in_dir)

    test_ht = khmer.load_counting_hash(os.path.join(in_dir, 'backup.ht'))
    test_good_read = 'CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'
    test_good_read2 = 'TAGTATCATCAAGGTTCAAGATGTTAATGAATAACAATTGCGCAGCAA'
    assert test_ht.count(test_good_read[:17]) > 0
    assert test_ht.count(test_good_read2[:17]) > 0
    
    assert status == 0
    assert os.path.exists(os.path.join(in_dir, 'backup.ht'))
    assert out.count('Backup: Saving') == 2
    assert 'Nothing' in out

def test_normalize_by_median_empty():
    CUTOFF='1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-empty.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

def test_normalize_by_median_empty_fq():
    CUTOFF='1'

    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-empty.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

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

def test_load_graph_no_tags():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = runscript(script, args)
    assert status == 0

    ht_file = outfile + '.ht'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    ht = khmer.load_hashbits(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...

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

def _DEBUG_make_graph(infilename, SIZE=1e7, N=2, K=20,
                do_partition=False,
                annotate_partitions=False,
                stop_big_traverse=False):
    script = scriptpath('load-graph.py')
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data(infilename)

    args.extend([outfile, infile])

    (status, out, err) = DEBUG_runscript(script, args)
    assert status == 0

    ht_file = outfile + '.ht'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    if do_partition:
	print ">>>> DEBUG: Partitioning <<<"
        script = scriptpath('partition-graph.py')
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        (status, out, err) = DEBUG_runscript(script, args)
        print out
        print err
        assert status == 0

	print ">>>> DEBUG: Merging Partitions <<<"
        script = scriptpath('merge-partitions.py')
        args = [outfile, '-k', str(K)]
        (status, out, err) = DEBUG_runscript(script, args)
        print out
        print err
        assert status == 0

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
	    print ">>>> DEBUG: Annotating Partitions <<<"
            script = scriptpath('annotate-partitions.py')
            args = ["-k", str(K), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            (status, out, err) = DEBUG_runscript(script, args, in_dir)
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

    print (status, out, err)

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '1001 2 98 1.0', line

def test_abundance_dist_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17, BIGCOUNT=False)

    script = scriptpath('abundance-dist.py')
    args = ['-z', htfile, infile, outfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    print (status, out, err)

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '255 2 98 1.0', line

def test_abundance_dist_single():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('abundance-dist-single.py')
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', infile, outfile]
    (status, out, err) = runscript(script, args, in_dir)
    print status
    print out
    print err
    assert status == 0

    print open(outfile).read()

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '1001 2 98 1.0', line

def test_abundance_dist_single_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('abundance-dist-single.py')
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-b', infile, outfile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    print (status, out, err)
    print open(outfile).read()

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '255 2 98 1.0', line

def test_do_partition():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('do-partition.py')
    args = ["-k", "20", graphbase, seqfile]

    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0, (out,err)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [ r.name.split('\t')[1] for r in screed.open(partfile) ]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1

def test_do_partition_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('do-partition.py')
    args = ["-k", "21", graphbase, seqfile]

    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0, (out,err)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [ r.name.split('\t')[1] for r in screed.open(partfile) ]
    parts = set(parts)

    assert len(parts) == 99, len(parts)

