DATA_DIR = '../data'
SCRIPTS_DIR = '../scripts'

###

import os
import glob
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

scriptsdir = os.path.join(thisdir, SCRIPTS_DIR)
scriptsdir = os.path.abspath(scriptsdir)

datadir = os.path.join(thisdir, DATA_DIR)
datadir = os.path.abspath(datadir)

import sys
import khmer
import subprocess


def test_quick_do_partition_calc():
    script = os.path.join(scriptsdir, 'do-partition-calc.py')
    datafile = os.path.join(datadir, '25k.fa')

    print 'running', script, datafile

    p = subprocess.Popen([sys.executable, script, datafile],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    print '---\n', out
    print '---\n', err
    assert p.returncode == 0


def test_quick_do_th_subset_calc():
    script = os.path.join(scriptsdir, 'do-th-subset-calc.py')
    datafile = os.path.join(datadir, '25k.fa')

    print 'running', script, datafile

    p = subprocess.Popen([sys.executable, script, datafile],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    print '---\n', out
    print '---\n', err
    assert p.returncode == 0


def test_quick_do_th_subset_calc():
    script = os.path.join(scriptsdir, 'do-th-subset-calc.py')
    datafile = os.path.join(datadir, '25k.fa')

    print 'running', script, datafile

    p = subprocess.Popen([sys.executable, script, datafile],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    print '---\n', out
    print '---\n', err
    assert p.returncode == 0


def test_quick_do_th_subset_save():
    script = os.path.join(scriptsdir, 'do-th-subset-save.py')
    datafile = os.path.join(datadir, '25k.fa')

    print 'running', script, datafile

    p = subprocess.Popen([sys.executable, script, datafile],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    print '---\n', out
    print '---\n', err
    assert p.returncode == 0


def test_quick_do_th_subset_load():
    script = os.path.join(scriptsdir, 'do-th-subset-load.py')
    datafile = os.path.join(datadir, '25k.fa')
    datapath = os.path.join(datadir, '25k.fa.*.pmap')
    files = glob.glob(datapath)

    x = [sys.executable, script, datafile]
    x.extend(files)

    print 'running', x

    p = subprocess.Popen(x, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    print '---\n', out
    print '---\n', err
    assert p.returncode == 0


def test_quick_graph_size_py():
    script = os.path.join(scriptsdir, 'graph-size-py.py')
    datafile = os.path.join(datadir, "*.gz")

    x = [sys.executable, script, 'occ.out', datafile]
    print 'running', x

    p = subprocess.Popen(x, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    print '---\n', out
    print '---\n', err
    assert p.returncode == 0
