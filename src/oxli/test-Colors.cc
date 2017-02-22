/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2013-2015, Michigan State University.
Copyright (C) 2015, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"
#include "oxli/hashbits.hh"
#include "oxli/labelhash.hh"
#include <iostream>

using namespace oxli;

int main()
{
    HashIntoType sizes[] = { 100000003, 100000004, 100000007, 10000000011};
    std::vector<HashIntoType> sizes_vec (sizes,
                                         sizes + sizeof(sizes) / sizeof(HashIntoType) );

    oxli::LabelHash * lh_pointer = new oxli::LabelHash(20, sizes_vec);
    oxli::Nodegraph * hb_pointer = (oxli::Hashbits *)lh_pointer;

    std::cout << "lh_pointer n_tags: " << lh_pointer->n_tags() << std::endl;
    std::cout << "hb_pointer n_tags: " << hb_pointer->n_tags() << std::endl;

    return 0;
}
