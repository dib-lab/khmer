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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Note that this file originally did not compile.  Me (holtgrew) change it
// it it compiled.  However, it might not make sense to fix it so it actually
// works.
// ==========================================================================

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_MISC_MAP_H
#define SEQAN_HEADER_MISC_MAP_H

#include <algorithm>
#include "base.h"


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    //

    template <typename TPair>
    inline typename TPair::T1 & keyOf(TPair &pair) {
        return getValueI1(pair);
    }
    template <typename TPair>
    inline typename TPair::T1 const & keyOf(TPair const &pair) {
        return getValueI1(pair);
    }
    template <typename TPair>
    inline typename TPair::T2 & objectOf(TPair &pair) {
        return getValueI2(pair);
    }
    template <typename TPair>
    inline typename TPair::T2 const & objectOf(TPair const &pair) {
        return getValueI2(pair);
    }

    //////////////////////////////////////////////////////////////////////////////
    //

//    template <typename TKey>
//    struct Map {
//        typedef std::map<TKey> Type;
//    };
    template <typename TKey, typename TObject>
    class Map< Pair<TKey, TObject> > {
//        typedef std::set< Pair<TKey, TObject>, SetLess< Pair<TKey, TObject> > > Type;
    };

}

#endif

