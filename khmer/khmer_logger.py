# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
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
"""Lightweight logging framework for khmer."""

import sys

__QUIET__ = False


def configure_logging(quiet):
    """Set the logging level."""
    global __QUIET__  # pylint: disable=global-statement
    __QUIET__ = quiet


def log_info(message, **kwargs):
    """For non-critical informative/status output to stderr."""
    if not __QUIET__:
        if kwargs:
            message = message.format(**kwargs)
        print(message, file=sys.stderr)


def log_error(message, **kwargs):
    """For critical error output to stderr."""
    if kwargs:
        message = message.format(**kwargs)
    print(message, file=sys.stderr)


def log_debug(message, **kwargs):
    """For non-critical debug output to stderr."""
    if not __QUIET__:
        if kwargs:
            message = message.format(**kwargs)
        print(message, file=sys.stderr)


def log_warn(message, **kwargs):
    """For critical warning output to stderr."""
    if kwargs:
        message = message.format(**kwargs)
    print(message, file=sys.stderr)


# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
