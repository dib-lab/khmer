#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
python ctb-iterative-bench-2.py
python filter-exact.py foo.fa xxx
python filter-inexact-all.py foo.fa yyy
python filter-inexact-any.py foo.fa zzz
python filter-inexact-run.py foo.fa aaa
python occupy.py foo.fa 12
