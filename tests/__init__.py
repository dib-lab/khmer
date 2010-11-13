import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

pythondir = os.path.join(thisdir, '..', 'python')
pythondir = os.path.abspath(pythondir)

import sys
sys.path.insert(0, pythondir)
