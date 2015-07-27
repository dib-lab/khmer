from __future__ import print_function, unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# Lightweight logging framework for khmer

import sys


class Logger(object):

    """
    Object to handle error/warning/info output in scripts.

    Four levels: info, warn, err, debug.

    Info is for status output, written to stderr, is surpressable.

    Err is for critical failures, will not be surpressed when quiet is set; Err
    output is written to stderr.

    Debug output is written to stdout and is surpressable.

    Warning output is written to stdout and is not surpressable.
    """

    def __init__(self, quiet):
        self._quiet = quiet

    def log_inf(self, message, **kwargs):
        """For non-critical informative/status output to stderr."""
        if not self._quiet:
            print(message.format(**kwargs), file=sys.stderr)

    def log_err(self, message, **kwargs):
        """For critical error output to stderr."""
        print(message.format(**kwargs), file=sys.stderr)

    def log_dbg(self, message, **kwagrs):
        """For non-critical debug output to stdout."""
        if not self._quiet:
            print(message.format(**kwargs), file=sys.stdout)

    def log_wrn(self, message, **kwargs):
        """For critical warning output to stdout."""
        print(message.format(**kwargs), file=sys.stdout)

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
