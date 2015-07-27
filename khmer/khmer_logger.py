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
    Object to handle error/warning/info output in scripts

    Four levels: info, warn, err, debug

    Info is for status output, written to stderr, is surpressable

    Warn is for warning output, written to stderr, is surpressable

    Err is for critical failures, will not be surpressed when quiet is set; Err
    output is written to stderr

    Debug output is written to stdout and is never surpressed
    """

    def __init__(self, quiet):
        self._quiet = quiet

    def log_inf(self, message, **kwargs):
        """For non-critical informative/status output"""
        if not self._quiet:
            print(message.format(**kwargs), file=sys.stderr)

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
