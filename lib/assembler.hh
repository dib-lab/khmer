/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2015-2016, The Regents of the University of California.

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
#ifndef ASSEMBLER_HH
#define ASSEMBLER_HH

#include <queue>
#include <functional>

#include "khmer.hh"

#include "khmer_exception.hh"
#include "read_parsers.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"


namespace khmer
{

class Hashtable;

// A function which takes a Kmer and returns true if it
// is to be filtered / ignored
typedef std::function<bool (Kmer&)> KmerFilter;
typedef std::list<KmerFilter> KmerFilterList;
// list instead of vector because: better insertation at ends, more
// space efficient, don't need random access
typedef std::list<Kmer> OrderedKmers;


inline bool apply_kmer_filters(Kmer& node, std::list<KmerFilter>& filters)
{
    if (!filters.size()) {
        return false;
    }

    for(auto filter : filters) {
        if (filter(node)) {
            return true;
        }
    }
    return false;
}


class AssemblerTraverser: public Traverser
{

protected:

    OrderedKmers contig_kmers;
    Kmer cursor;
    KmerFilterList filters;
    const bool traverse_right; // 0 for left, 1 for right

public:

    explicit AssemblerTraverser(const Hashtable * ht,
                             Kmer start_kmer,
                             KmerFilterList filters,
                             bool traverse_right = 1);

    void gather_linear_path();
    unsigned int get_path_length() const;
    std::string build_contig() const;
    std::string assemble();


};

class Assembler: public Traverser
{
    friend class Hashtable;

public:

    explicit Assembler(const Hashtable * ht);

    std::string assemble_linear_path(const Kmer seed_kmer,
                                     const Hashtable * stop_bf=0) const;

    std::string _assemble_right(const std::string start_kmer,
                                std::list<KmerFilter>& node_filters) const;

    std::string _assemble_left(const std::string start_kmer,
                               std::list<KmerFilter>& node_filters) const;

    //std::string _assemble_directed(const char * start_kmer,
    //                               const Hashtable * stop_bf,
    //                               const bool assemble_left=false) const;

};

}
#endif
