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

#include <functional>

#include "khmer.hh"
#include "kmer_hash.hh"

#include "labelhash.hh"
#include "traversal.hh"
#include "kmer_filters.hh"


namespace khmer
{

class Hashtable;
class LabelHash;

/**
 * \class LinearAssembler
 *
 * \brief Naively assemble linear paths.
 *
 * The LinearAssember is the basic building block of general assembly.
 * The core function, _assemble_directed, takes a direction and branch-specialiazed
 * AssemblerTraverser and gathers a contig up until any branch points. The behavior
 * of this function can be modified by passing it a cursor with varying
 * filtering functions.
 *
 * The assemble, assemble_right, and assemble_left functions are for convenience.
 * They take care of building the filters and creating the contig strings.
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
class LinearAssembler
{
    friend class Hashtable;
    const Hashtable * graph;
    WordLength _ksize;

public:

    explicit LinearAssembler(const Hashtable * ht);

    std::string assemble(const Kmer seed_kmer,
                         const Hashtable * stop_bf = 0) const;

    std::string assemble_right(const Kmer seed_kmer,
                        const Hashtable * stop_bf = 0) const;

    std::string assemble_left(const Kmer seed_kmer,
                        const Hashtable * stop_bf = 0) const;

    template <bool direction>
    std::string _assemble_directed(AssemblerTraverser<direction>& cursor) const;
};

// The explicit specializations need to be declared in the same translation unit
// as their unspecialized declaration.
template<>
std::string LinearAssembler::_assemble_directed<LEFT>(AssemblerTraverser<LEFT>& cursor) const;

template<>
std::string LinearAssembler::_assemble_directed<RIGHT>(AssemblerTraverser<RIGHT>& cursor) const;


class LabeledLinearAssembler
{
    friend class Hashtable;
    const LinearAssembler * linear_asm;
    const Hashtable * graph;
    const LabelHash * lh;
    WordLength _ksize;

public:

    explicit LabeledLinearAssembler(const LabelHash * lh);

    StringVector assemble(const Kmer seed_kmer,
                         const Hashtable * stop_bf=0) const;

    template <bool direction>
    void _assemble_directed(NonLoopingAT<direction>& start_cursor,
                            StringVector& paths) const;

};

}
#endif
