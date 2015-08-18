from __future__ import print_function
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
import tempfile
import os
import shutil
from pkg_resources import Requirement, resource_filename, ResolutionError
import nose
import sys
import traceback
import subprocess
from io import open


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


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

cleanup_list = []


def get_temp_filename(filename, tempdir=None):
    if tempdir is None:
        tempdir = tempfile.mkdtemp(prefix='khmertest_')
        cleanup_list.append(tempdir)

    return os.path.join(tempdir, filename)


def cleanup():
    global cleanup_list

    for path in cleanup_list:
        shutil.rmtree(path, ignore_errors=True)
    cleanup_list = []


def scriptpath(scriptname='interleave-reads.py'):
    "Return the path to the scripts, in both dev and install situations."

    # note - it doesn't matter what the scriptname is here, as long as
    # it's some khmer script present in this version of khmer.

    path = os.path.join(os.path.dirname(__file__), "../scripts")

    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname, sandbox=False):
    """
    Find & run a script with exec (i.e. not via os.system or subprocess).
    """

    import pkg_resources
    ns = {"__name__": "__main__"}
    ns['sys'] = globals()['sys']

    try:
        pkg_resources.get_distribution("khmer").run_script(
            scriptname, ns)
        return 0
    except pkg_resources.ResolutionError as err:
        if sandbox:
            path = os.path.join(os.path.dirname(__file__), "../sandbox")
        else:
            path = scriptpath()

        scriptfile = os.path.join(path, scriptname)
        if os.path.isfile(scriptfile):
            if os.path.isfile(scriptfile):
                exec(compile(open(scriptfile).read(), scriptfile, 'exec'), ns)
                return 0
        elif sandbox:
            raise nose.SkipTest("sandbox tests are only run in a repository.")

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
        except nose.SkipTest:
            raise
        except SystemExit as e:
            status = e.code
        except:
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
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()

        out = out.decode('utf-8')
        err = err.decode('utf-8')

        if p.returncode != 0 and not fail_ok:
            print('out:', out)
            print('err:', err)
            raise AssertionError("exit code is non zero: %d" % p.returncode)

        return (p.returncode, out, err)
    finally:
        os.chdir(cwd)


def longify(listofints):
    """List of ints => list of longs, only on py2.

    Takes a list of numeric types, and returns longs on python2, or the
    original list on python3.
    """
    # For map(long, [list of ints]) cross-version hackery
    if sys.version_info.major < 3:
        return map(long, listofints)
    return listofints
