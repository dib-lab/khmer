#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import glob

filelist = glob.glob('*R1*.fastq.gz')

for r1 in filelist:
    r2 = r1.replace('R1', 'R2')
    final_pe = r1[:-9] + '.pe.fq.gz'
    final_se = r1[:-9] + '.se.fq.gz'
    print("""\
mkdir trim
cd trim
java -jar /usr/local/bin/trimmomatic-0.30.jar PE ../%s ../%s s1_pe s1_se s2_pe s2_se ILLUMINACLIP:/usr/local/share/adapters/TruSeq3-PE.fa:2:30:10
/usr/local/share/khmer/scripts/interleave-reads.py s1_pe s2_pe | gzip -9c > ../%s

cat s1_se s2_se | gzip -9c > ../%s
cd ..
rm -r ./trim/

chmod u-w %s %s
""" % (r1, r2, final_pe, final_se, final_pe, final_se))
