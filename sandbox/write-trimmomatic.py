#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import glob
import sys

filelist = glob.glob('*R1*.fastq.gz')

for r1 in filelist:
    r2 = r1.replace('R1', 'R2')
    final_pe = r1[:-9] + '.pe.fq.gz'
    final_se = r1[:-9] + '.se.fq.gz'
    print """\
mkdir trim
cd trim
TrimmomaticPE ../%s ../%s s1_pe s1_se s2_pe s2_se ILLUMINACLIP:%s:2:30:10
interleave-reads.py s1_pe s2_pe | gzip -9c > ../%s

cat s1_se s2_se | gzip -9c > ../../trimmed/%s
cd ..
rm -r ./trim/

chmod u-w %s %s
""" % (r1, r2, sys.argv[1], final_pe, final_se, final_pe, final_se)
