#! /usr/bin/env python
import sys, time
import subprocess

# velveth %(filename)s.ass $K -short -paired %(filename)s

#COMMAND="velvetg %(filename)s.ass -read_trkg yes -exp_cov 3 -cov_cutoff 0 -min_contig_lgth 1000"
COMMAND="/root/khmer/velvet-assemble.sh %(filename)s 31"


N_PROCESSES = 8

filenames = sys.argv[1:]

running = []

while 1:
    while len(running) < N_PROCESSES:
        filename = filenames.pop(0)
        cmd = COMMAND % dict(filename=filename)
        p = subprocess.Popen(cmd, shell=True)
        print 'running:', cmd
        running.append(p)

    i = 0;
    while i < len(running):
        p = running[i]
        returncode = p.poll()
        if returncode is not None:
            print 'done! removing', p
            running.remove(p)
        else:
            i += 1

    time.sleep(1)

print '\n\n** done **\n'
