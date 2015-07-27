from __future__ import print_function, unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# Lightweight logging framework for khmer

import sys

global __QUIET__
__QUIET__ = False


def configure_logging(quiet):
    global __QUIET__
    __QUIET__ = quiet


def log_info(message, **kwargs):
    """For non-critical informative/status output to stderr."""
    global __QUIET__
    if not __QUIET__:
        if kwargs:
            message = message.format(**kwargs)
        print(message, file=sys.stderr)


def log_error(message, **kwargs):
    """For critical error output to stderr."""
    if kwargs:
        message = message.format(**kwargs)
    print(message, file=sys.stderr)


def log_debug(message, **kwagrs):
    """For non-critical debug output to stderr."""
    global __QUIET__
    if not __QUIET__:
        if kwargs:
            message = message.format(**kwargs)
        print(message, file=sys.stderr)


def log_warn(message, **kwargs):
    """For critical warning output to stderr."""
    if kwargs:
        message = message.format(**kwargs)
    print(message, file=sys.stderr)


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
