from __future__ import print_function
import sys
sys.path.insert(0, '../')
import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = '../khmer/_version.py'
versioneer.versionfile_build = '../khmer/_version.py'
versioneer.tag_prefix = 'v'  # tags are like v1.2.0
versioneer.parentdir_prefix = '..'

print(versioneer.get_version())
