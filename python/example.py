#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
sys.path.insert(0, 'build/lib.linux-i686-2.3/')

import khmer
ktable = khmer.new_ktable(6)
ktable.consume("ATGAGAGACACAGGGAGAGACCCAATTAGAGAATTGGACC")
for i in range(0, ktable.n_entries()):
    n = ktable.get(i)
    if n:
        print ktable.reverse_hash(i), "is present", n, "time(s)."
