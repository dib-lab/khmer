#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import os
import sys
import time
import subprocess

localdir = os.path.dirname(__file__)
localdir = os.path.abspath(localdir)

COMMAND = "%(localdir)s/stats-assemble.sh %(localdir)s %(filename)s"


N_PROCESSES = 8

print '** running with %d concurrent processes' % N_PROCESSES

filenames = sys.argv[1:]
filenames.reverse()

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

            print 'done! %d of ~%d' % (report_number + 1, report_number +
                                       len(running) + len(filenames))
            report_number += 1
            running.remove((cmd, p))
        else:
            i += 1

    time.sleep(1)

print '\n** done **\n'

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
