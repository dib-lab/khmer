# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
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
# pylint: disable=C0111,C0103
"""
Test code that verifies current script output md5 hashes against recorded
hashes, to ensure that script function isn't changing.
"""
from __future__ import print_function
from __future__ import absolute_import

import khmer
from . import khmer_tst_utils as utils

#
# hashes recorded as of git commit 799039ffcf15d2a3ac6902ae62ae2da81030e8d2
# (for trim-low-abund without --diginorm behavior) and
# b939a34b565ce973224abdd0eeb53d6b52833c01 (for trim-low-abund updates).
#


def test_normalize_by_median_k21_C20_M1e7():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('normalize-by-median.py', ['-C', '20', '-k', '21',
                                               '-M', '1e7', '-o', outfile,
                                               infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == '5fed8ed10896639ad7fb607f93a4473d', hash


def test_normalize_by_median_k21_C15_M1e7():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('normalize-by-median.py', ['-C', '15', '-k', '21',
                                               '-M', '1e7', '-o', outfile,
                                               infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == 'fd38b9e17a736694c20e66a3d8411eaa', hash


def test_trim_low_abund_k21_C0_M1e7_diginorm():
    # should be same as normalize-by-median -C 20 -k 21 -M 1e7
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '0', '-k', '21',
                                          '--diginorm',
                                          '--diginorm-coverage', '20',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == '5fed8ed10896639ad7fb607f93a4473d', hash


def test_trim_low_abund_k21_C0_M1e7_diginorm_dn15():
    # should be same as normalize-by-median -C 15 -k 21 -M 1e7
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '0', '-k', '21',
                                          '--diginorm',
                                          '--diginorm-coverage', '15',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == 'fd38b9e17a736694c20e66a3d8411eaa', hash


def test_trim_low_abund_k21_C2_M1e7_diginorm_dn15():
    # should be same as normalize-by-median -C 15 -k 21 -M 1e7
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '2', '-k', '21',
                                          '--diginorm',
                                          '--diginorm-coverage', '15',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == 'e57be1402a9fd84d8481d8ef5a695596', hash


def test_trim_low_abund_k21_M1e7_C2():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '2', '-k', '21',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == '21ad2f64ceeb1e153f1e513486df9434', hash


def test_trim_low_abund_k21_M1e7_C3():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '3', '-k', '21',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == '3da84a42d8538e24240e78cd0eadf168', hash


def test_trim_low_abund_k21_M1e7_C4():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '4', '-k', '21',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == '578f124937b5a639985cc0a808139653', hash


def test_trim_low_abund_k21_M1e7_C4_variable():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '4', '-k', '21',
                                          '-V',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == '154d2a25599b3dcd2cf9cd2375999e56', hash


def test_trim_low_abund_k21_M1e7_C4_variable_Z25():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '4', '-k', '21',
                                          '-V', '-Z', '25',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    hash = utils._calc_md5(open(outfile, 'rb'))
    assert hash == 'eb8bc6ed6246fb819100b9e7dde6fad2', hash


def test_trim_low_abund_k21_M1e7_C4_variable_Z15():
    infile = utils.get_test_data('simple-genome-reads.fa')
    outfile = utils.get_temp_filename('out')
    utils.runscript('trim-low-abund.py', ['-C', '4', '-k', '21',
                                          '-V', '-Z', '15',
                                          '-M', '1e7', '-o', outfile,
                                          infile])

    with open(outfile, 'rb') as output:
        hashval = utils._calc_md5(output)
    assert hashval == 'b948ae4f2114be27279de37b46376f6d', hash
