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

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_H

#include <seqan/index/pizzachili_api.h>
#include <seqan/index/index_pizzachili_string.h>

namespace SEQAN_NAMESPACE_MAIN {

/*!
 * @defgroup PizzaChiliIndexFibres Pizza &amp; Chili Index Fibres
 *
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        PizzaChiliIndex @endlink index.
 *
 * Pizza &amp; Chili indices are compressed indices. Hence, this fibre is used for
 * searching in the index.
 *
 * @see Fibre
 * @see Index#getFibre
 * @see PizzaChiliIndex
 *
 * @tag PizzaChiliIndexFibres#PizzaChiliText
 * @brief The original text the index is based on.
 *
 * @tag PizzaChiliIndexFibres#PizzaChiliCompressed
 * @brief The compressed suffix array.
 */

struct FibrePizzaChiliCompressed_;

typedef Tag<FibreText_> const FibrePizzaChiliText;
typedef Tag<FibrePizzaChiliCompressed_> const FibrePizzaChiliCompressed;

typedef FibrePizzaChiliText PizzaChiliText;
typedef FibrePizzaChiliCompressed PizzaChiliCompressed;

/*!
 * @class PizzaChiliIndex Pizza & Chili Index
 *
 * @extends Index
 *
 * @headerfile <seqan/index.h>
 *
 * @brief An adapter for the Pizza &amp; Chili index API.
 *
 * @signature template <typename TText, typename TSpec>
 *            class Index<TText, PizzaChili<TSpec> >;
 *
 * @tparam TSpec Tag specifying the Pizza &amp; Chili index library to use. Types:
 *               @link PizzaChiliIndexTags @endlink
 * @tparam TText The text type. Types: @link String @endlink
 *
 * @see PizzaChiliString
 * @see PizzaChiliIndexFibres
 * @see IndexFindAlgorithm
 * @see PizzaChiliIndexTags
 */

template <typename TText, typename TSpec>
class Index<TText, PizzaChili<TSpec> > {
public:
    typedef typename Value<TText>::Type TValue;
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

    impl::index_t index_handle;
    Holder<String<TValue, PizzaChili<TSpec> > > text;

    Index() : index_handle(0), text() { }

    Index(Index& other) : index_handle(0), text() {
SEQAN_CHECKPOINT
        // Explicitly request the other's index text.
        setIndexText(*this, indexText(other));
    }

    Index(Index const& other) : index_handle(0), text() {
SEQAN_CHECKPOINT
        // Explicitly request the other's index text.
        setIndexText(*this, indexText(other));
    }

    template <typename TOtherText>
    Index(TOtherText& txt) : index_handle(0), text() {
SEQAN_CHECKPOINT
        setIndexText(*this, txt);
    }

    ~Index() {
SEQAN_CHECKPOINT
        clear(*this);
    }

    Index& operator =(Index const& other) {
SEQAN_CHECKPOINT
        if (this == &other)
            return *this;

        clear(*this);
        setIndexText(*this, indexText(other));

        return *this;
    }

private:
    //Index(Index const& other);
    //Index operator =(Index const& other);
};

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChiliText> {
    typedef String<typename Value<TText>::Type, PizzaChili<TSpec> >& Type;
};

template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChiliText> {
    typedef String<typename Value<TText>::Type, PizzaChili<TSpec> > const& Type;
};

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TText, typename TSpec>
    inline void
    clearIndex(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

        if (me.index_handle != 0) {
            impl::error_t e =
                TCodeProvider::free_index(me.index_handle);

            if (e != 0)
                SEQAN_ABORT(TCodeProvider::error_index(e));

            me.index_handle = 0;
        }
    }
} // namespace impl

template <typename TText, typename TSpec>
inline void
clear(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    impl::clearIndex(me);
    clear(me.text);
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {

    template <typename TText, typename TSpec>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<TSpec> >& /*me*/) {
SEQAN_CHECKPOINT
        return "";
    }

    template <typename TText>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<PizzaChiliSA> >& /*me*/) {
SEQAN_CHECKPOINT
        return "copy_text";
    }

    template <typename TText>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<PizzaChiliFM> >& /*me*/) {
SEQAN_CHECKPOINT
        return "-a 0";
    }

    template <typename TText>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<PizzaChiliRsa> >& /*me*/) {
SEQAN_CHECKPOINT
        return "copy_text";
    }

} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChiliText>::Type
indexText(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    return getFibre(me, PizzaChiliText());
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChiliText>::Type
indexText(Index<TText, PizzaChili<TSpec> > const& me) {
SEQAN_CHECKPOINT
    return getFibre(me, PizzaChiliText());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChiliText>::Type
getFibre(Index<TText, PizzaChili<TSpec> > const& me, PizzaChiliText const) {
SEQAN_CHECKPOINT
    return value(me.text);
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChiliText>::Type
getFibre(Index<TText, PizzaChili<TSpec> >& me, PizzaChiliText const) {
SEQAN_CHECKPOINT
    return value(me.text);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool
indexSupplied(Index<TText, PizzaChili<TSpec> >& me, PizzaChiliCompressed const) {
SEQAN_CHECKPOINT
    return me.index_handle != 0;
}

template <typename TText, typename TSpec>
inline bool
indexSupplied(Index<TText, PizzaChili<TSpec> >& me, PizzaChiliText const) {
SEQAN_CHECKPOINT
    return length(value(me.text)) > 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool
indexSolveDependencies(Index<TText, PizzaChili<TSpec> >& me, PizzaChiliCompressed const) {
SEQAN_CHECKPOINT
    return indexSupplied(me, PizzaChiliText());
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TText, typename TSpec>
    inline bool
    createPizzaChiliIndex(
        Index<TText, PizzaChili<TSpec> >& me,
        uchar_t* textstart,
        ulong_t textlength
    ) {
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
        // Read-only access, therefore safe cast.
        char* options = const_cast<char*>(impl::getOptionsString(me));
        impl::error_t e =
            TCodeProvider::build_index(textstart, textlength, options, &me.index_handle);

        if (e != 0) {
            // TODO(holtgrew): Do we need a logging interface?
            //SEQAN_REPORT(TCodeProvider::error_index(e));
            //SEQAN_REPORT(options);
            me.index_handle = 0;
            return false;
        }

        value(me.text) = String<typename Value<TText>::Type, PizzaChili<TSpec> >(me.index_handle);

        return true;
    }
} // namespace impl

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, PizzaChili<TSpec> >& me, PizzaChiliCompressed const) {
    typedef
        typename RemoveConst_<
            typename Index<TText, PizzaChili<TSpec> >::TValue
        >::Type alph_t;

    SEQAN_ASSERT_EQ(sizeof(alph_t), 1);
    SEQAN_ASSERT((IsSameType<typename IsSimple<alph_t>::Type, True>::VALUE));

    impl::clearIndex(me);

    impl::uchar_t* textstart =
        reinterpret_cast<impl::uchar_t*>(
            const_cast<alph_t*>(indexText(me).data_begin)
        );
    impl::ulong_t textlength = length(indexText(me));
    return impl::createPizzaChiliIndex(me, textstart, textlength);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec, typename TOtherText>
inline void
setIndexText(Index<TText, PizzaChili<TSpec> >& me, TOtherText& text) {
    //typedef
    //    typename RemoveConst_<
    //        typename Value<TOtherText>::Type
    //    >::Type alph_t;
    clear(me);

    /*
    SEQAN_ASSERT(IsContiguous<TOtherText>::VALUE)
    SEQAN_ASSERT_EQ(BitsPerValue<alph_t>::VALUE, 8)
    SEQAN_ASSERT((IsSameType<typename IsSimple<alph_t>::Type, True>::VALUE));

    String<alph_t, CStyle> cstr = text;
    impl::uchar_t* textstart =
        reinterpret_cast<impl::uchar_t*>(
            const_cast<alph_t*>(static_cast<alph_t const*>(cstr))
        );
    impl::ulong_t textlength = length(text);
    impl::createPizzaChiliIndex(me, textstart, textlength);
    */
    getFibre(me, PizzaChiliText()) = text;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool open(
    Index<TText, PizzaChili<TSpec> >& me,
    char const* filename
) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
    clear(me);
    impl::error_t e =
        TCodeProvider::load_index(const_cast<char*>(filename), &me.index_handle);
    if (e != 0) {
        SEQAN_REPORT(TCodeProvider::error_index(e));
    }
    else
        value(me.text) = String<typename Value<TText>::Type, PizzaChili<TSpec> >(me.index_handle);
    return e == 0;
}

template <typename TText, typename TSpec>
inline bool save(
    Index<TText, PizzaChili<TSpec> >& me,
    char const* filename
) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
    // Before we can save the index we have to construct it, in case it isn't
    // already constructed.
    indexRequire(me, PizzaChiliCompressed());
    impl::error_t e =
        TCodeProvider::save_index(me.index_handle, const_cast<char*>(filename));
    if (e != 0) {
        SEQAN_REPORT(TCodeProvider::error_index(e));
    }
    return e == 0;
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_H
