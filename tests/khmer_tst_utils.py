# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
# Copyright (C) 2016, Google, Inc
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
# pylint: disable=missing-docstring

import tempfile
import os
import shutil
import pkg_resources
from pkg_resources import Requirement, resource_filename, ResolutionError
import sys
import traceback
import subprocess
from io import open  # pylint: disable=redefined-builtin
from hashlib import md5

from khmer import reverse_complement as revcomp

import pytest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def _equals_rc(query, match):
    return (query == match) or (revcomp(query) == match)


def _contains_rc(match, query):
    return (query in match) or (revcomp(query) in match)


def _calc_md5(fp):
    m = md5()
    m.update(fp.read())
    return m.hexdigest()


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(
            Requirement.parse("khmer"), "khmer/tests/test-data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                                filename)
    return filepath


CLEANUPLIST = []


def get_temp_filename(filename, tempdir=None):
    if tempdir is None:
        tempdir = tempfile.mkdtemp(prefix='khmertest_')
    CLEANUPLIST.append(tempdir)

    return os.path.join(tempdir, filename)


def cleanup():
    global CLEANUPLIST  # pylint: disable=global-statement

    for path in CLEANUPLIST:
        shutil.rmtree(path, ignore_errors=True)
    CLEANUPLIST = []


def scriptpath(scriptname='interleave-reads.py'):
    """Return the path to the scripts, in both dev and install situations."""
    # note - it doesn't matter what the scriptname is here, as long as
    # it's some khmer script present in this version of khmer.

    path = os.path.join(os.path.dirname(__file__), "../scripts")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    path = os.path.join(os.path.dirname(__file__), "../../EGG-INFO/scripts")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname, sandbox=False):
    """Find & run a script with exec (i.e. not via os.system or subprocess)."""
    namespace = {"__name__": "__main__"}
    namespace['sys'] = globals()['sys']

    try:
        pkg_resources.get_distribution("khmer").run_script(
            scriptname, namespace)
        return 0
    except pkg_resources.ResolutionError:
        pass

    if sandbox:
        path = os.path.join(os.path.dirname(__file__), "../sandbox")
    else:
        path = scriptpath()

    scriptfile = os.path.join(path, scriptname)
    if os.path.isfile(scriptfile):
        if os.path.isfile(scriptfile):
            exec(compile(open(scriptfile).read(), scriptfile, 'exec'),
                 namespace)
            return 0
    else:
        raise RuntimeError("Tried to execute {} but it is"
                           " not a file.".format(scriptfile))

    return -1


def runscript(scriptname, args, in_directory=None,
              fail_ok=False, sandbox=False):
    """Run a Python script using exec().

    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.

    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    sysargs = [scriptname]
    sysargs.extend(args)

    cwd = os.getcwd()

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        if in_directory:
            os.chdir(in_directory)
        else:
            in_directory = cwd

        try:
            print('running:', scriptname, 'in:', in_directory, file=oldout)
            print('arguments', sysargs, file=oldout)

            status = _runscript(scriptname, sandbox=sandbox)
        except SystemExit as err:
            status = err.code
        except:  # pylint: disable=bare-except
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        os.chdir(cwd)

    if status != 0 and not fail_ok:
        print(out)
        print(err)
        assert False, (status, out, err)

    return status, out, err


def run_shell_cmd(cmd, fail_ok=False, in_directory=None):
    cwd = os.getcwd()
    if in_directory:
        os.chdir(in_directory)

    print('running: ', cmd)
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        (out, err) = proc.communicate()

        out = out.decode('utf-8')
        err = err.decode('utf-8')

        if proc.returncode != 0 and not fail_ok:
            print('out:', out)
            print('err:', err)
            raise AssertionError("exit code is non zero: %d" % proc.returncode)

        return (proc.returncode, out, err)
    finally:
        os.chdir(cwd)


def longify(listofints):
    """List of ints => list of longs, only on py2.

    Takes a list of numeric types, and returns longs on python2, or the
    original list on python3.
    """
    # For map(long, [list of ints]) cross-version hackery
    if sys.version_info.major < 3:
        return map(long, listofints)  # pylint: disable=bad-builtin
    return listofints


def copy_test_data(testfile, newfilename=None):
    basename = os.path.basename(testfile)
    if newfilename is not None:
        basename = newfilename
    infile = get_temp_filename(basename)
    shutil.copyfile(get_test_data(testfile), infile)
    return infile
