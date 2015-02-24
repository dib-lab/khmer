// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Konrad Ludwig Moritz Rudolph <konrad.rudolph@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_FIND_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_FIND_H

namespace SEQAN_NAMESPACE_MAIN {

struct PizzaChiliFinder_;

typedef Tag<PizzaChiliFinder_> const PizzaChiliFinder;

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
struct DefaultFinder<Index<TText, PizzaChili<TSpec> > > {
    typedef PizzaChiliFinder Type;
};

//////////////////////////////////////////////////////////////////////////////

// We define a special suffix array fibre for our finder. This fibre is used to
// represent the search results and our search results aren't actually part of
// a suffix array but rather a C array of ulong values.
template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> >, Tag<FibreSA_> const> {
    typedef impl::ulong_t* Type;
};

//////////////////////////////////////////////////////////////////////////////

// The following specialization of the finder is necessary only because we have
// to free the data from the C API in the destructor.
template <typename TText, typename TSpec>
class Finder<Index<TText, PizzaChili<TSpec> >, PizzaChiliFinder>
    : public Finder<Index<TText, PizzaChili<TSpec> >, Default>
{
    typedef Finder<Index<TText, PizzaChili<TSpec> >, Default> TBase;

    typedef typename TBase::TIndex TIndex;
    typedef typename TBase::TSA TSA;
    typedef typename TBase::TIterator TIterator;

public:

    using TBase::index;
    using TBase::range;
    using TBase::data_iterator;

    Finder() : TBase() { }

    Finder(TIndex& index) : TBase(index) { }

    Finder(TIndex const& index) : TBase(index) { }

    ~Finder() {
SEQAN_CHECKPOINT
        if (range.i1 != 0)
            std::free(range.i1);
    }
};

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TPattern>
    inline uchar_t* getPizzaChiliString(TPattern const& pattern) {
SEQAN_CHECKPOINT
        typedef
            typename RemoveConst_<
                typename Value<TPattern>::Type
            >::Type alph_t;

        // This const_cast is safe because 'cstr' is only read.
        // On the other hand, this cast is necessary to prevent a copy from
        // being made: this would result in an invalid pointer to local memory.
        String<alph_t, CStyle> const cstr = *const_cast<TPattern*>(&pattern);
        // This const_cast is safe because the return value is only read.
        return
            reinterpret_cast<uchar_t*>(
                const_cast<alph_t*>(static_cast<alph_t const*>(cstr))
            );
    }

    inline uchar_t* getPizzaChiliString(char const* pattern) {
SEQAN_CHECKPOINT
        return reinterpret_cast<uchar_t*>(const_cast<char*>(pattern));
    }
} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec, typename TSpecFinder, typename TPattern>
inline void _findFirstIndex(
    Finder<Index<TText, PizzaChili<TSpec> >, TSpecFinder>& finder,
    TPattern const& pattern,
    PizzaChiliFinder const
) {
SEQAN_CHECKPOINT
    typedef Index<TText, PizzaChili<TSpec> > TIndex;
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

    if (finder.range.i1 != 0)
        std::free(finder.range.i1);

    TIndex& index = haystack(finder);
    indexRequire(index, PizzaChiliCompressed());

    impl::uchar_t* c_pattern =
        impl::getPizzaChiliString(pattern);
    impl::ulong_t patternlength = length(pattern);

    impl::ulong_t numocc;
    impl::ulong_t* occ;

    impl::error_t e =
        TCodeProvider::locate(
            index.index_handle,
            c_pattern,
            patternlength,
            &occ,
            &numocc
        );

    if (e != 0)
        SEQAN_ABORT(TCodeProvider::error_index(e));

    finder.range.i1 = occ;
    finder.range.i2 = occ + numocc;
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_FIND_H
