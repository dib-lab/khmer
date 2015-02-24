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

#ifndef SEQAN_HEADER_FIND_MULTI_H
#define SEQAN_HEADER_FIND_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct MultiPatternFinder_;
typedef Tag<MultiPatternFinder_> MultipatternFinder;

//____________________________________________________________________________

template <typename THaystack>
class Finder<THaystack, MultipatternFinder>
{
//____________________________________________________________________________
private:
    unsigned int data_pattern;

public:
    Finder():
        data_pattern(0)
    {
SEQAN_CHECKPOINT
    }

    Finder(Finder const & other_):
        data_pattern(other_.data_pattern)
    {
SEQAN_CHECKPOINT
    }

    ~Finder()
    {
SEQAN_CHECKPOINT
    }
//____________________________________________________________________________

    Finder &
    operator = (Finder const & other_)
    {
SEQAN_CHECKPOINT
        data_pattern = other_.data_pattern;
        return *this;
    }
//____________________________________________________________________________



//////////////////////////////////////////////////////////////////////////////
/*
template <typename THaystack, typename TNeedle>
bool
findNext(Finder & me,
         THaystack & hstk,
         TNeedle const & ndl)
{
SEQAN_CHECKPOINT
    ++hstk;
    return find(me, hstk, ndl);
}*/

//////////////////////////////////////////////////////////////////////////////

};

template <typename THaystack>
inline unsigned int &
needle(Finder<THaystack, MultipatternFinder> & me)
{
SEQAN_CHECKPOINT
    return me.data_pattern;
}
template <typename THaystack>
inline unsigned int const &
needle(Finder<THaystack, MultipatternFinder> const & me)
{
SEQAN_CHECKPOINT
    return me.data_pattern;
}
//____________________________________________________________________________

template <typename THaystack>
inline void
setNeedle(Finder<THaystack, MultipatternFinder> & me, unsigned int const needleIndex_)
{
SEQAN_CHECKPOINT
    me.data_pattern = needleIndex_;
}

//____________________________________________________________________________

template <typename THaystack>
inline void
init(Finder<THaystack, MultipatternFinder> & me)
{
SEQAN_CHECKPOINT
    me.data_pattern = 0;
}
//____________________________________________________________________________


template <typename THaystack, typename THaystack2, typename TNeedle>
inline bool
find(Finder<THaystack, MultipatternFinder> & me,
     THaystack2 & hstk,
     TNeedle const & ndl)
{
SEQAN_CHECKPOINT
    while ( needle(me) < length(ndl) )
    {
        Finder<THaystack2, Horspool> horspool(ndl[needle(me)]);
        bool found = find(horspool, hstk, ndl[needle(me)]);
        if (found)
        {
            return true;
        }
        setPosition(hstk, 0);
        ++needle(me);
    }
    return false;
}
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
