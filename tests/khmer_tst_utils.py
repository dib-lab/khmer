import tempfile, os, shutil

thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

tempfile.tempdir = thisdir

def get_test_data(filename):
    return os.path.join(thisdir, 'test-data', filename)

cleanup_list = []

def get_temp_filename(filename, tempdir=None):
    if tempdir is None:
        tempdir = tempfile.mkdtemp(prefix='khmertest_')
        cleanup_list.append(tempdir)
    
    return os.path.join(tempdir, filename)

def cleanup():
    global cleanup_list

    for path in cleanup_list:
        shutil.rmtree(path)
    cleanup_list = []
