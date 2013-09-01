#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import tempfile
from screed.fasta import fasta_iter
import shutil
import os.path
import subprocess


def load_sequences(filename):
    d = {}
    records = list(fasta_iter(open(filename), parse_description=False))

    for r in records:
        name = r['name']
        partition = name.rsplit('\t', 1)[1]
        partition = int(partition)

        x = d.get(partition, [])
        x.append(r)
        d[partition] = x

    return len(records), d


def assemble_sequences(records, k, length_cutoff=1000):
    dirname = tempfile.mkdtemp()
    os.chdir(dirname)

    try:
        seqfile = os.path.join(dirname, 'seqs.fa')
        fp = open(seqfile, 'w')
        for r in records:
            fp.write('>%s\n%s\n' % (r['name'].split()[0], r['sequence']))
        fp.close()

        p = subprocess.Popen(
            ['python', '/root/khmer/scripts/strip-and-split-for-assembly.py',
             'seqs.fa seqs.fa'], shell=True)
        p.communicate()
        assert p.returncode == 0

        assemble_dir = os.path.join(dirname, 'assemble')
        p = subprocess.Popen(
            'velveth %s %d -shortPaired %s.pe -short %s.se' % (
            assemble_dir, k, seqfile, seqfile),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        p = subprocess.Popen(
            'velvetg %s -read_trkg yes -exp_cov auto -cov_cutoff 0' % (
                assemble_dir,),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        x = []
        total = 0
        for r in fasta_iter(open(os.path.join(assemble_dir, 'contigs.fa'))):
            seqlen = len(r['sequence'])
            if seqlen >= length_cutoff:
                x.append(r)
                total += seqlen

        return total, x
    finally:
        shutil.rmtree(dirname)
        # print 'XXX', dirname


def best_assemble_sequences(r, try_k=(33, 35, 37, 39, 41, 43, 45, 47, 49, 51)):

    best_k = try_k[0]
    best_total, best_records = assemble_sequences(r, best_k)

    for k in try_k[1:]:
        total, records = assemble_sequences(r, k)
        if total > best_total:
            best_total = total
            best_records = records
            best_k = k

    return best_k, best_total, best_records

n, partitions = load_sequences(sys.argv[1])

print 'loaded %d sequences in %d partitions' % (n, len(partitions))

fp = open(os.path.basename(sys.argv[1]) + '.best', 'w')
for pid in partitions:
    print 'trying', pid
    records = partitions[pid]
    k, total, records = best_assemble_sequences(records)
    print 'best assembly for part %d: k=%d, %d bp' % (pid, k, total)

    for n, r in enumerate(records):
        fp.write('>part%d.%d\n%s\n' % (pid, n, r['sequence']))

fp.close()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
