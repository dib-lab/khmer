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
#ifndef GMAP_HH
#define GMAP_HH

#include <functional>
#include <memory>

#include "oxli.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashgraph.hh"

namespace oxli {


template <typename T, bool threadsafe>
class GuardedHashMap {

    public:

        // Filter should be owned exclusively by GuardedKmerMap
        std::unique_ptr<Nodegraph> filter;
        std::map<HashIntoType, T> data;
        
        explicit GuardedHashMap(WordLength ksize,
                                unsigned short n_tables,
                                uint64_t max_table_size)
        {
            std::vector<uint64_t> table_sizes = get_n_primes_near_x(n_tables, max_table_size);
            filter = std::unique_ptr<Nodegraph>(new Nodegraph(ksize, table_sizes));
        }

        T get(HashIntoType kmer) const
        {
            if (filter->get_count(kmer)) {
                auto search = data.find(kmer);
                if (search != data.end()) {
                    return search->second;
                }
            }
            
            return NULL;
        }

        void set(HashIntoType kmer, T item)
        {
            filter->count(kmer);
            data[kmer] = item;
        }

        bool contains(HashIntoType kmer) const 
        {
            return get(kmer) != NULL;
        }

        uint64_t size() const 
        {
            return data.size();
        }
};

template <typename T>
class GuardedHashMap<T, true>: public GuardedHashMap<T, false>
{
    private:

        uint32_t lock;

    public:

        using GuardedHashMap<T, false>::GuardedHashMap;
        using GuardedHashMap<T, false>::filter;
        using GuardedHashMap<T, false>::data;

        explicit GuardedHashMap(WordLength ksize,
                                unsigned short n_tables,
                                uint64_t max_table_size) : 
            GuardedHashMap<T, false>(ksize, n_tables, max_table_size),
            lock(0)
        {
        }
        
        T get(HashIntoType kmer) const 
        {
            if (filter->get_count(kmer)) {
                while(!__sync_bool_compare_and_swap( &lock, 0, 1));
                auto search = data.find(kmer);
                if (search != data.end()) {
                    __sync_bool_compare_and_swap( &lock, 1, 0);
                    return search->second;
                }
                __sync_bool_compare_and_swap( &lock, 1, 0);
            }

            return NULL;
        }

        void set(HashIntoType kmer, T item) 
        {
            while(!__sync_bool_compare_and_swap( &lock, 0, 1));
            set(kmer, item);
            __sync_bool_compare_and_swap( &lock, 1, 0);
        }
};

}

#endif
