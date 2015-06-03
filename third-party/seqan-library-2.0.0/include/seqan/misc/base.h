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

#ifndef SEQAN_HEADER_MISC_BASE_H
#define SEQAN_HEADER_MISC_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // In SeqAn sets and maps store elements as pairs of (key,object)
    // the elements of sets without objects are the keys.
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // Key/Object meta-functions
    //////////////////////////////////////////////////////////////////////////////

/* moved to basic_h, see #6
    template <typename TElement>
    struct Key {
        typedef TElement Type;
    };
*/
/* moved to basic_aggregates
    template <typename TKey, typename TObject, typename TSpec>
    struct Key< Pair<TKey, TObject, TSpec> >
    {
        typedef TKey Type;
    };
*/

//////////////////////////////////////////////////////////////////////////////

    template <typename TElement>
    struct Object {
        typedef Nothing Type;
    };

    template <typename TKey, typename TObject, typename TSpec>
    struct Object< Pair<TKey, TObject, TSpec> > {
        typedef TObject Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // keyOf function
    //////////////////////////////////////////////////////////////////////////////

    template <typename TElement>
    inline typename Key<TElement const>::Type &
    keyOf(TElement const & element)
    {
        return element;
    }

    template <typename TKey, typename TObject, typename TSpec>
    inline TKey const &
    keyOf(Pair<TKey, TObject, TSpec> const &element) {
        return element.i1;
    }

    //////////////////////////////////////////////////////////////////////////////
    // objectOf
    //////////////////////////////////////////////////////////////////////////////

    template <typename TKey, typename TObject, typename TSpec>
    inline TObject &
    objectOf(Pair<TKey, TObject, TSpec> &element) {
        return element.i2;
    }
    template <typename TKey, typename TObject, typename TSpec>
    inline TObject const &
    objectOf(Pair<TKey, TObject, TSpec> const &element) {
        return element.i2;
    }

}

#endif

