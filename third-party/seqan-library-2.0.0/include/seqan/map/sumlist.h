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

#ifndef SEQAN_HEADER_SUMLIST_H
#define SEQAN_HEADER_SUMLIST_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//Default Specialization Tag

template <typename TSpec = Default>
struct SkipSumList;

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec = SkipSumList<> >
class SumList;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct DIMENSION
{
    enum { VALUE = 1}; //dummy implementation to make VC++ happy
};

template <unsigned int DIM_, typename TValue, typename TSpec>
struct DIMENSION< SumList<DIM_, TValue, TSpec> >
{
    enum { VALUE = DIM_};
};
template <unsigned int DIM_, typename TValue, typename TSpec>
struct DIMENSION< SumList<DIM_, TValue, TSpec> const>
{
    enum { VALUE = DIM_};
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
struct Value< SumList<DIM, TValue, TSpec> >
{
    typedef TValue Type;
};
template <unsigned int DIM, typename TValue, typename TSpec>
struct Value< SumList<DIM, TValue, TSpec> const >
{
    typedef TValue Type;
};


//////////////////////////////////////////////////////////////////////////////
// SumListValues:
// Stores one DIM tupel of Values
//////////////////////////////////////////////////////////////////////////////


template <unsigned int DIM, typename TValue>
struct SumListValues
{
    TValue values[DIM];

    SumListValues(MinimalCtor)
    {}
    SumListValues()
    {
        clear(*this);
    }
    SumListValues(SumListValues const & other)
    {
        *this = other;
    }
    SumListValues(TValue const * arr)
    {
        *this = arr;
    }
    ~SumListValues()
    {
    }
    SumListValues const & operator = (SumListValues const & other)
    {
        arrayCopy(other.values, other.values + DIM, values);
        return *this;
    }
    SumListValues const & operator = (TValue const * arr)
    {
        arrayCopy(arr, arr + DIM, values);
        return *this;
    }

    TValue & operator [] (unsigned int pos)
    {
        return values[pos];
    }
    TValue const & operator [] (unsigned int pos) const
    {
        return values[pos];
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Values;

template <unsigned int DIM, typename TValue, typename TSpec>
struct Values< SumList<DIM, TValue, TSpec> >
{
    typedef SumListValues<DIM, TValue> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline SumListValues<DIM, TValue> const &
operator += (SumListValues<DIM, TValue> & left,
             SumListValues<DIM, TValue> const & right)
{
    for (unsigned int i = 0; i < DIM; ++i)
    {
        left.values[i] += right.values[i];
    }
    return left;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline SumListValues<DIM, TValue> const &
operator -= (SumListValues<DIM, TValue> & left,
             SumListValues<DIM, TValue> const & right)
{
    for (unsigned int i = 0; i < DIM; ++i)
    {
        left.values[i] -= right.values[i];
    }
    return left;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline bool
operator == (SumListValues<DIM, TValue> const & left,
             SumListValues<DIM, TValue> const & right)
{
    for (unsigned int i = 0; i < DIM; ++i)
    {
        if (left.values[i] != right.values[i]) return false;
    }
    return true;
}
template <unsigned int DIM, typename TValue>
inline bool
operator != (SumListValues<DIM, TValue> const & left,
             SumListValues<DIM, TValue> const & right)
{
    for (unsigned int i = 0; i < DIM; ++i)
    {
        if (left.values[i] != right.values[i]) return true;
    }
    return false;
}
//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline void
clear(SumListValues<DIM, TValue> & me)
{
    arrayFill(me.values, me.values + DIM, 0);
}


//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
