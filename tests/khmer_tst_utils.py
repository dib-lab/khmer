#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import tempfile
import os
import shutil
from pkg_resources import Requirement, resource_filename, ResolutionError
from cStringIO import StringIO
import nose
import sys
import khmer.file


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


def _runscript(scriptname, sandbox):
    import pkg_resources
    ns = {"__name__": "__main__"}
    ns['sys'] = globals()['sys']

    try:
        pkg_resources.get_distribution("khmer").run_script(
            scriptname, ns)
        return 0
    except pkg_resources.ResolutionError, err:
        if sandbox:
            paths = [os.path.join(os.path.dirname(__file__), "../sandbox")]
        else:
            paths = [os.path.join(os.path.dirname(__file__),
                                  "../scripts")]
            paths.extend(os.environ['PATH'].split(':'))
        for path in paths:
            scriptfile = os.path.join(path, scriptname)
            if os.path.isfile(scriptfile):
                execfile(scriptfile, ns)
                return 0
        if sandbox:
            raise nose.SkipTest("sandbox tests are only run in a repository.")

    return -1


def runscript(scriptname, args, in_directory=None, fail_ok=False, sandbox=False):
    """
    Run the given Python script, with the given args, in the given directory,
    using 'execfile'.
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
        sys.stderr = StringIO()

        if in_directory:
            os.chdir(in_directory)

        try:
            status = _runscript(scriptname, sandbox)
        except nose.SkipTest:
            raise
        except SystemExit, e:
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
        print out
        print err
        assert False, (status, out, err)

    return status, out, err
