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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_QUASAR_H
#define SEQAN_HEADER_FIND_QUASAR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Quasar
//////////////////////////////////////////////////////////////////////////////

//DISABLED .Class.Pattern.param.TSpec.type:Spec.Quasar

struct Quasar_;
typedef Tag<Quasar_> Quasar;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, Quasar> {
//____________________________________________________________________________
private:
    Pattern(Pattern const& other);
    Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
    Holder<TNeedle> data_host;


//____________________________________________________________________________

    Pattern() {
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
SEQAN_CHECKPOINT
        setHost(*this, ndl);
    }

    ~Pattern() {
        SEQAN_CHECKPOINT
    }
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, Quasar> >
{
    typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, Quasar> const>
{
    typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
inline void
setHost (Pattern<TNeedle, Quasar> & me, TNeedle2 const& needle)
{
    SEQAN_CHECKPOINT
    setValue(me.data_host, needle);
}

template <typename TNeedle, typename TNeedle2>
inline void
setHost (Pattern<TNeedle, Quasar> & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, Quasar> & /*me*/)
{
SEQAN_CHECKPOINT
}


//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Quasar>const>::Type &
host(Pattern<TNeedle, Quasar> & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Quasar>const>::Type &
host(Pattern<TNeedle, Quasar> const & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool
find(TFinder & finder, Pattern<TNeedle, Quasar> & me)
{
    SEQAN_CHECKPOINT

    if (empty(finder)) {
        _patternInit(me);
        _finderSetNonEmpty(finder);
    } else {
        finder+=4;
    }

    return true;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H

