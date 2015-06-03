// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains Index specializations for thrust::device_vector.
// ==========================================================================

#ifndef SEQAN_INDEX_DEVICE_H_
#define SEQAN_INDEX_DEVICE_H_

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Device                                                  [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct Device<Index<TText, TSpec> >
{
    typedef Index<typename Device<TText>::Type, TSpec>          Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                                [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct IsDevice<Index<TText, TSpec> > : IsDevice<TText> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                                 [Iter]
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct IsDevice<Iter<TIndex, VSTree<TSpec> > > : IsDevice<TIndex> {};

// ----------------------------------------------------------------------------
// Metafunction FibreSA                                          [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type> Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreLcp                                         [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreLcp>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreLcp>
{
    typedef thrust::device_vector<typename Size<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreChildtab                                    [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreChildtab>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreChildtab>
{
    typedef thrust::device_vector<typename Size<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreBwt                                         [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>   Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSA                                        [Device FMIndex]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec, typename TConfig>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TSpec, TConfig> >, FibreSA>
{
    typedef CompressedSA<thrust::device_vector<TValue, TAlloc>, TSpec, TConfig>      Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec, typename TConfig>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, FMIndex<TSpec, TConfig> >, FibreSA>
{
    typedef CompressedSA<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec, TConfig>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibrePrefixSums                                     [Device LF]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec, typename TConfig>
struct Fibre<LF<thrust::device_vector<TValue, TAlloc>, TSpec, TConfig>, FibrePrefixSums>
{
    typedef thrust::device_vector<TValue, TAlloc>               TText_;
    typedef typename Size<LF<TText_, TSpec, TConfig> >::Type    TSize_;
//    typedef typename Value<LF<TText_, TSpec, TConfig> >::Type  TValue_;
//    typedef Tuple<TSize_, ValueSize<TValue_>::VALUE>                Type;

    typedef thrust::device_vector<TSize_>                       Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec, typename TConfig>
struct Fibre<LF<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec, TConfig>, FibrePrefixSums>
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec> TText_;
    typedef typename Size<LF<TText_, TSpec, TConfig> >::Type            TSize_;
//    typedef typename Value<LF<TText_, TSpec, TConfig> >::Type           TValue_;
//    typedef Tuple<TSize_, ValueSize<TValue_>::VALUE>                    Type;

    typedef thrust::device_vector<TSize_>                               Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreBwt                                            [Device LF]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec, typename TConfig>
struct Fibre<LF<thrust::device_vector<TValue, TAlloc>, TSpec, TConfig>, FibreBwt>
{
    typedef thrust::device_vector<TValue, TAlloc>               TText_;
    typedef typename Value<LF<TText_, TSpec, TConfig> >::Type   TValue_;

    typedef RankDictionary<TValue_, TwoLevels<Device<TSpec> > > Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec, typename TConfig>
struct Fibre<LF<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec, TConfig>, FibreBwt>
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec> TText_;
    typedef typename Value<LF<TText_, TSpec, TConfig> >::Type           TValue_;

    typedef RankDictionary<TValue_, TwoLevels<Device<TSpec> > >         Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSentinels                                      [Device LF]
// ----------------------------------------------------------------------------

// NOTE(esiragusa): Single text sentinel rank dictionary is an integer.
//template <typename TText, typename TViewSpec, typename TSpec>
//struct Fibre<LF<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSentinels>
//{
//    typedef typename Fibre<LF<TText, TSpec, TConfig>, FibreSentinels>::Type   Type;
//};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec, typename TConfig>
struct Fibre<LF<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec, TConfig>, FibreSentinels>
{
//    typedef RankDictionary<bool, TwoLevels<Device<TSpec> > >    Type;
    typedef RankDictionary<bool, Naive<Device<TSpec> > >    Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreRanks                              [Device RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TValue, TwoLevels<Device<TSpec> > >, FibreRanks>
{
    typedef thrust::device_vector<RankDictionaryEntry_<TValue, TwoLevels<TSpec> > > Type;
};

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TValue, Naive<Device<TSpec> > >, FibreRanks>
{
    typedef RankDictionary<TValue, Naive<TSpec> >               TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type               TSize_;

    typedef thrust::device_vector<TSize_>                       Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSparseString                         [Device CompressedSA]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec, typename TConfig>
struct Fibre<CompressedSA<thrust::device_vector<TValue, TAlloc>, TSpec, TConfig>, FibreSparseString>
{
    typedef thrust::device_vector<TValue/*, TAlloc*/>       TText_;
    typedef typename SAValue<TText_>::Type                  TSAValue_;
    typedef thrust::device_vector<TSAValue_/*, TAlloc*/>    TString_;

    typedef SparseString<TString_, TSpec>                   Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec, typename TConfig>
struct Fibre<CompressedSA<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec, TConfig>, FibreSparseString>
{
    typedef StringSet<thrust::device_vector<TValue/*, TAlloc*/>, TSSetSpec> TText_;
    typedef typename SAValue<TText_>::Type                                  TSAValue_;
    typedef thrust::device_vector<TSAValue_/*, TAlloc*/>                    TString_;

    typedef SparseString<TString_, TSpec>                                   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreIndicators                           [Device SparseString]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<SparseString<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreIndicators>
{
    typedef RankDictionary<bool, TwoLevels<Device<TSpec> > >        Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreValues                               [Device SparseString]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): FibreValues is already defined as below.

//template <typename TValue, typename TAlloc, typename TSpec>
//struct Fibre<SparseString<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreValues>
//{
//    typedef thrust::device_vector<TValue, TAlloc>       Type;
//};

// ----------------------------------------------------------------------------
// Metafunction HistoryStack_                              [Device VSTree Iter]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>, VSTree<TSpec> > >
{
    typedef thrust::device_vector<TValue, TAlloc>       TText_;
    typedef Index<TText_, TIndexSpec>                   TIndex_;
    typedef Iter<TIndex_, VSTree<TSpec> >               TIter_;
    typedef typename HistoryStackEntry_<TIter_>::Type   TEntry_;
    typedef thrust::device_vector<TEntry_>              Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>, VSTree<TSpec> > >
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>     TText_;
    typedef Index<TText_, TIndexSpec>                                       TIndex_;
    typedef Iter<TIndex_, VSTree<TSpec> >                                   TIter_;
    typedef typename HistoryStackEntry_<TIter_>::Type                       TEntry_;
    typedef thrust::device_vector<TEntry_>                                  Type;
};

// ============================================================================
// Functions
// ============================================================================
// NOTE(esiragusa): Functions assign() are not specific to Device Index.

// ----------------------------------------------------------------------------
// Function assign()                                                  [FMIndex]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TText2, typename TOccSpec2, typename TSpec2>
inline void
assign(Index<TText, FMIndex<TSpec, TConfig> > & index, Index<TText2, FMIndex<TOccSpec2, TSpec2> > const & source)
{
    assign(indexText(index), indexText(source));
    assign(indexLF(index), indexLF(source));
    assign(indexSA(index), indexSA(source));

    // Set the pointer.
    setFibre(indexSA(index), indexLF(index), FibreLF());
}

template <typename TText, typename TSpec, typename TConfig, typename TText2, typename TOccSpec2, typename TSpec2>
inline void
assign(Index<TText, FMIndex<TSpec, TConfig> > & index, Index<TText2, FMIndex<TOccSpec2, TSpec2> > & source)
{
    assign(index, reinterpret_cast<Index<TText2, FMIndex<TOccSpec2, TSpec2> > const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                                       [LF]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TText2, typename TSpec2, typename TConfig>
inline void
assign(LF<TText, TSpec, TConfig> & lf, LF<TText2, TSpec2, TConfig> const & source)
{
    assign(getFibre(lf, FibrePrefixSums()), getFibre(source, FibrePrefixSums()));
    assign(getFibre(lf, FibreBwt()), getFibre(source, FibreBwt()));
    assign(getFibre(lf, FibreSentinels()), getFibre(source, FibreSentinels()));
    assign(lf.sentinelSubstitute, source.sentinelSubstitute);
}

template <typename TText, typename TSpec, typename TText2, typename TSpec2, typename TConfig>
inline void
assign(LF<TText, TSpec, TConfig> & lf, LF<TText2, TSpec2, TConfig> & source)
{
    assign(lf, reinterpret_cast<LF<TText2, TSpec2, TConfig> const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TValue, TwoLevels<TSpec> > & dict, RankDictionary<TValue2, TwoLevels<TSpec2> > const & source)
{
    assign(getFibre(dict, FibreRanks()), getFibre(source, FibreRanks()));
    assign(dict._length, source._length);
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TValue, TwoLevels<TSpec> > & dict, RankDictionary<TValue2, TwoLevels<TSpec2> > & source)
{
    assign(dict, reinterpret_cast<RankDictionary<TValue2, TwoLevels<TSpec2> > const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TValue, Naive<TSpec> > & dict, RankDictionary<TValue2, Naive<TSpec2> > const & source)
{
    assign(getFibre(dict, FibreRanks()), getFibre(source, FibreRanks()));
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TValue, Naive<TSpec> > & dict, RankDictionary<TValue2, Naive<TSpec2> > & source)
{
    assign(dict, reinterpret_cast<RankDictionary<TValue2, Naive<TSpec2> > const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                             [CompressedSA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TText2, typename TSpec2, typename TConfig>
inline void
assign(CompressedSA<TText, TSpec, TConfig> & sa, CompressedSA<TText2, TSpec2, TConfig> const & source)
{
    assign(getFibre(sa, FibreSparseString()), getFibre(source, FibreSparseString()));
}

template <typename TText, typename TSpec, typename TText2, typename TSpec2, typename TConfig>
inline void
assign(CompressedSA<TText, TSpec, TConfig> & sa, CompressedSA<TText2, TSpec2, TConfig> & source)
{
    assign(getFibre(sa, FibreSparseString()), getFibre(source, FibreSparseString()));
}

// ----------------------------------------------------------------------------
// Function assign()                                             [SparseString]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
inline void
assign(SparseString<TString, TSpec> & sparseString, SparseString<TString2, TSpec2> const & source)
{
    assign(getFibre(sparseString, FibreValues()), getFibre(source, FibreValues()));
    assign(getFibre(sparseString, FibreIndicators()), getFibre(source, FibreIndicators()));
    assign(sparseString._length, source._length);
}

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
inline void
assign(SparseString<TString, TSpec> & sparseString, SparseString<TString2, TSpec2> & source)
{
    assign(sparseString, reinterpret_cast<SparseString<TString2, TSpec2> const &>(source));
}

}

#endif  // #ifndef SEQAN_INDEX_DEVICE_H_
