import sys
sys.path.insert(0, 'build/lib.linux-i686-2.3/')

import khmer
ktable = khmer.new_ktable(6)
ktable.consume("ATGAGAGACACAGGGAGAGACCCAATTAGAGAATTGGACC")
for i in range(0, ktable.n_entries()):
    n = ktable.get(i)
    if n:
        print ktable.reverse_hash(i), "is present", n, "time(s)."
