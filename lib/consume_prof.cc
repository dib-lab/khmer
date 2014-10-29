//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#include "khmer.hh"
#include "hashtable.hh"

int main(int argc, char *argv[])
{
    unsigned int total_reads;
    unsigned long long n_consumed;

    khmer::Hashtable ht(20, 16777216 + 1);
    ht.consume_fasta(argv[1], total_reads, n_consumed);
}
