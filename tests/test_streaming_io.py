from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from . import khmer_tst_utils as utils
import subprocess
import os.path
import difflib


def scriptpath():
    path = os.path.join(os.path.dirname(__file__), "../scripts")
    return path


def files_are_equal(a, b):
    al = open(a).readlines()
    bl = open(b).readlines()

    return al == bl


def diff_files(a, b):
    al = open(a).readlines()
    bl = open(b).readlines()

    results = "\n".join(difflib.context_diff(al, bl, fromfile=a, tofile=b))
    return results


def run_shell_cmd(cmd, fail_ok=False):
    print('running: ', cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    (out, err) = p.communicate()

    if p.returncode != 0 and not fail_ok:
        raise AssertionError("exit code is non zero: %d" % p.returncode)

    return (p.returncode, out, err)


def test_interleave_split_1():
    in1 = utils.get_test_data('paired.fq.1')
    in2 = utils.get_test_data('paired.fq.2')
    
    out1 = utils.get_temp_filename('a.fa')
    out2 = utils.get_temp_filename('b.fa')

    cmd = """
       {scripts}/interleave-reads.py {in1} {in2} -o -             |
       {scripts}/split-paired-reads.py -1 {out1} -2 {out2}
    """

    cmd = cmd.format(scripts=scriptpath(),
                     in1=in1, in2=in2,
                     out1=out1, out2=out2)

    run_shell_cmd(cmd)

    assert files_are_equal(in1, out1), diff_files(in1, out1)
    assert files_are_equal(in2, out2), diff_files(in2, out2)


def test_extract_paired_pe():
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.pe')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py -p - -s /dev/null > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_extract_paired_se():
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.se')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py -p /dev/null -s - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_norm_by_median_1():
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.pe')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py -p - -s /dev/null |
       {scripts}/normalize-by-median.py -p - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)
