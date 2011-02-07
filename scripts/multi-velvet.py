#! /usr/bin/env python
import os
import sys, time
import subprocess

K=33

localdir = os.path.dirname(__file__)
localdir = os.path.abspath(localdir)

# velveth %(filename)s.ass $K -short -paired %(filename)s

#COMMAND="velvetg %(filename)s.ass -read_trkg yes -exp_cov 3 -cov_cutoff 0 -min_contig_lgth 1000"
COMMAND="%%(localdir)s/velvet-assemble.sh %%(localdir)s %%(filename)s %d" % K


N_PROCESSES = 8
print '** K is %d' % K
print '** running with %d concurrent processes' % N_PROCESSES

filenames = sys.argv[1:]

running = []

report_number = 0
while filenames or running:
    while filenames and len(running) < N_PROCESSES:
        filename = filenames.pop(0)
        cmd = COMMAND % dict(filename=filename, localdir=localdir)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        print 'running:', cmd
        running.append((cmd, p))

    i = 0
    while i < len(running):
        (cmd, p) = running[i]
        returncode = p.poll()
        if returncode is not None:
            out = p.stdout.read()
            err = p.stderr.read()

            report_fp = open('multi-velvet.report.%d' % report_number, 'w')
            report_fp.write('Command: %s' % cmd)
            report_fp.write('\nOutput:\n%s\n---' % out)
            report_fp.write('\nError:\n%s\n---' % err)

            print 'done! %d of ~%d' % (report_number + 1, report_number + len(running) + len(filenames))
            report_number += 1
            running.remove((cmd, p))
        else:
            i += 1

    time.sleep(1)

print '\n** done **\n'
