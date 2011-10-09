import sys, os, shutil
from cStringIO import StringIO
import traceback

import khmer_tst_utils as utils
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
    args = ['-C', '1', counting_ht, infile]
    (status, out, err) = runscript(script, args, in_dir)
    assert status == 0

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([ r.sequence for r in screed.open(outfile) ])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs
