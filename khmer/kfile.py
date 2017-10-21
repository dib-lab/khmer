# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014-2015, Michigan State University.
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
"""File handling/checking utilities for command-line scripts."""


import os
import sys
import errno
from stat import S_ISBLK, S_ISFIFO, S_ISCHR
import gzip
import bz2file


def check_input_files(file_path, force):
    """Check the status of the file.

    If the file is empty or doesn't exist AND if the file is NOT a
    fifo/block/named pipe then a warning is printed and sys.exit(1) is
    called
    """
    mode = None

    if file_path == '-':
        return

    try:
        mode = os.stat(file_path).st_mode
    except OSError:
        print("ERROR: Input file %s does not exist" %
              file_path, file=sys.stderr)

        if not force:
            print("NOTE: This can be overridden using the --force argument",
                  file=sys.stderr)
            print("Exiting", file=sys.stderr)
            sys.exit(1)
        else:
            return

    # block devices/stdin will be nonzero
    if S_ISBLK(mode) or S_ISFIFO(mode) or S_ISCHR(mode):
        return

    if not os.path.exists(file_path):
        print("ERROR: Input file %s does not exist; exiting" %
              file_path, file=sys.stderr)
        if not force:
            print("NOTE: This can be overridden using the --force argument",
                  file=sys.stderr)
            sys.exit(1)
    else:
        if os.stat(file_path).st_size == 0:
            print("ERROR: Input file %s is empty; exiting." %
                  file_path, file=sys.stderr)
            if not force:
                print("NOTE: This can be overridden using the --force"
                      " argument", file=sys.stderr)
                sys.exit(1)


def check_file_writable(file_path):
    """Return if file_path is writable, exit out if it's not."""
    try:
        file_obj = open(file_path, "a")
    except IOError as error:
        if error.errno == errno.EACCES:
            print("ERROR: File %s does not have write "
                  % file_path + "permission; exiting", file=sys.stderr)
            sys.exit(1)
        else:
            print("ERROR: " + error.strerror, file=sys.stderr)
    else:
        file_obj.close()
        return


def check_space(in_files, force, _testhook_free_space=None):
    """Check for available disk space.

    Estimate size of input files passed, then calculate disk space
    available and exit if disk space is insufficient.
    """
    # Get disk free space in Bytes assuming non superuser
    # and assuming all inFiles are in same disk
    in_file = in_files[0]

    dir_path = os.path.dirname(os.path.realpath(in_file))
    target = os.statvfs(dir_path)

    if _testhook_free_space is None:
        free_space = target.f_frsize * target.f_bavail
    else:
        free_space = _testhook_free_space

    # Check input file array, remove corrupt files
    valid_files = [f for f in in_files if os.path.isfile(f)]

    # Get input file size as worst case estimate of
    # output file size
    file_sizes = [os.stat(f).st_size for f in valid_files]
    total_size = sum(file_sizes)

    size_diff = total_size - free_space
    if size_diff > 0:
        mem_needed = '{:.1f} GB'.format(size_diff / 1e9)
        mem_requested = '{:.1f} GB'.format(total_size / 1e9)
        mem_available = '{:.1f} GB'.format(free_space / 1e9)
        message = 'Not enough free space on disk for output files;'
        message += '\n       Need at least {:s} more.'.format(mem_needed)
        message += '\n       Estimated Output size: {:s}'.format(mem_requested)
        message += '\n       Free space: {:s}'.format(mem_available)
        if force:
            print('WARNING:', message, file=sys.stderr)
        else:
            message = 'ERROR: ' + message + \
                '\nNOTE: This can be overridden using the --force argument'
            raise SystemExit(message)


def check_space_for_graph(outfile_name, hash_size, force,
                          _testhook_free_space=None):
    """
    Check that we have enough size to write the specified graph.

    The "hash_size" parameter should equal the total bytes required for the
    entire data structure.
    """
    dir_path = os.path.dirname(os.path.realpath(outfile_name))
    target = os.statvfs(dir_path)

    if _testhook_free_space is None:
        free_space = target.f_frsize * target.f_bavail
    else:
        free_space = _testhook_free_space  # allow us to test this code...

    size_diff = hash_size - free_space
    if size_diff > 0:
        mem_needed = '{:.1f} GB'.format(size_diff / 1e9)
        mem_requested = '{:.1f} GB'.format(hash_size / 1e9)
        mem_available = '{:.1f} GB'.format(free_space / 1e9)
        message = 'Not enough free space on disk for saved graph files;'
        message += '\n       Need at least {:s} more.'.format(mem_needed)
        message += '\n       Table size: {:s}'.format(mem_requested)
        message += '\n       Free space: {:s}'.format(mem_available)
        if force:
            print('WARNING:', message, file=sys.stderr)
        else:
            message = 'ERROR: ' + message + \
                '\nNOTE: This can be overridden using the --force argument'
            raise SystemExit(message)


def check_valid_file_exists(in_files):
    """Warn if input files are empty or missing.

    In a scenario where we expect multiple input files and
    are OK with some of them being empty or non-existent,
    this check warns to stderr if any input file is empty
    or non-existent.
    """
    for in_file in in_files:
        if in_file == '-':
            pass
        elif os.path.exists(in_file):
            mode = os.stat(in_file).st_mode
            if os.stat(in_file).st_size > 0 or S_ISBLK(mode) or S_ISFIFO(mode):
                return
            else:
                print('WARNING: Input file %s is empty' %
                      in_file, file=sys.stderr)
        else:
            print('WARNING: Input file %s not found' %
                  in_file, file=sys.stderr)


def is_block(fthing):
    """Take in a file object and checks to see if it's a block or fifo."""
    if fthing is sys.stdout or fthing is sys.stdin:
        return True
    else:
        mode = os.stat(fthing.name).st_mode
        return S_ISBLK(mode) or S_ISCHR(mode)


def describe_file_handle(fthing):
    """Return the name of file or a description."""
    if is_block(fthing):
        return "block device"
    else:
        return fthing.name


def add_output_compression_type(parser):
    """Add compression arguments to a parser object."""
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--gzip', default=False, action='store_true',
                       help='Compress output using gzip')
    group.add_argument('--bzip', default=False, action='store_true',
                       help='Compress output using bzip2')


def get_file_writer(file_handle, do_gzip, do_bzip):
    """Generate and return a file object with specified compression."""
    ofile = None

    if do_gzip and do_bzip:
        raise ValueError("Cannot specify both bzip and gzip compression!")

    if do_gzip:
        ofile = gzip.GzipFile(fileobj=file_handle, mode='w')
    elif do_bzip:
        ofile = bz2file.open(file_handle, mode='w')
    else:
        ofile = file_handle

    return ofile
