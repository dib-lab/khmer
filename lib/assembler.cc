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
#include "hashtable.hh"
#include "traversal.hh"
#include "assembler.hh"

using namespace khmer;
using namespace std;

Assembler::Assembler(const Hashtable * ht) :
    Traverser(ht)
{

}

// Starting from the given seed k-mer, assemble the maximal linear path in
// both directions.
//
// No guarantees on direction, of course - this may return the reverse
// complement of the input sequence.
//
// Note: as written, will ignore branches to the left and continue
// past them; this probably needs to be fixed.  For now, this means
// that assembling from two different directions may yield different
// results.
std::string Assembler::assemble_linear_path(const Kmer seed_kmer,
                                            const Hashtable * stop_bf)
    const
{
    std::string start_kmer = seed_kmer.get_string_rep(_ksize);
    std::string right = _assemble_right(start_kmer.c_str(), stop_bf);

    start_kmer = _revcomp(start_kmer);
    std::string left = _assemble_right(start_kmer.c_str(), stop_bf);

    left = left.substr(_ksize);
    return _revcomp(left) + right;
}


std::string Assembler::_assemble_right(const char * start_kmer,
                                       const Hashtable * stop_bf)
    const
{
    const char bases[] = "ACGT";
    std::string kmer = start_kmer;
    std::string contig = kmer;

    while (1) {
        const char * base = &bases[0];
        bool found = false;
        char found_base;
        bool found2 = false;

        while(*base != 0) {
            std::string try_kmer = kmer.substr(1) + (char) *base;

            // a hit!
            if (graph->get_count(try_kmer.c_str()) &&
                (!stop_bf || !stop_bf->get_count(try_kmer.c_str()))) {
                if (found) {
                    found2 = true;
                    break;
                }
                found_base = (char) *base;
                found = true;
            }
            base++;
        }
        if (!found || found2) {
            break;
        } else {
            contig += found_base;
            kmer = kmer.substr(1) + found_base;
            found = true;
        }
    }
    return contig;
}

/*
std::string Assembler::_assemble_directed(const char * start_kmer,
                                       const Hashtable * stop_bf,
                                       const bool assemble_left)
    const
{
    Kmer kmer = this->build_kmer(start_kmer);
    std::cout << "starting on kmer " << kmer.get_string_rep(_ksize) << std::endl;
    std::string contig = start_kmer;
    Traverser traverser(this);
    KmerQueue neighbors;
    unsigned short found;

    auto keep_func = [&] (Kmer& node) {
        return !stop_bf || !stop_bf->get_count(node);
    };

    std::cout << "start loop" << std::endl;
    while (1) {
        if (assemble_left) {
            found = traverser.traverse_left(kmer, neighbors, keep_func, 1);
        } else {
            found = traverser.traverse_right(kmer, neighbors, keep_func, 1);
        }
        //std::cout << "check traverser result" << std::endl;
        if (found == 0) {
            std::cout << "no neighbors, break" << std::endl;
            break;
        } else if (found > 1) {
            KmerQueue().swap(neighbors); // clear queue
            std::cout << "break" << std::endl;
            break;
        } else {
            //std::cout << "put base on contig" << std::endl;
            //contig += revtwobit_repr(kmer & 3);
            //std::cout << "get new kmer" << std::endl;
            kmer = neighbors.front();
            if (assemble_left) {
                contig += (kmer.get_string_rep(this->_ksize)[0]);
            } else {
                contig += (kmer.get_string_rep(this->_ksize)[this->_ksize-1]);

            }

            //std::cout << "pop!" << std::endl;
            neighbors.pop();
        }
    }
    return contig;
}
*/
