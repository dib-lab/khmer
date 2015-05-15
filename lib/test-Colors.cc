//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2013-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "khmer.hh"
#include "hashtable.hh"
#include "hashbits.hh"
#include "labelhash.hh"
#include <iostream>

using namespace khmer;

int main()
{
    HashIntoType sizes[] = { 100000003, 100000004, 100000007, 10000000011};
    std::vector<HashIntoType> sizes_vec (sizes,
                                         sizes + sizeof(sizes) / sizeof(HashIntoType) );

    khmer::LabelHash * lh_pointer = new khmer::LabelHash(20, sizes_vec);
    khmer::Hashbits * hb_pointer = (khmer::Hashbits *)lh_pointer;

    std::cout << "lh_pointer n_tags: " << lh_pointer->n_tags() << std::endl;
    std::cout << "hb_pointer n_tags: " << hb_pointer->n_tags() << std::endl;

    return 0;
}
