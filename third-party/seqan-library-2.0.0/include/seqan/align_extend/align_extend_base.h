// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// This file contains routines to extend an existing Align object
// ==========================================================================


#ifndef INCLUDE_ALIGN_ALIGN_EXTEND_BASE_H
#define INCLUDE_ALIGN_ALIGN_EXTEND_BASE_H

namespace seqan {


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = void>
struct AlignExtend_
{
};

// ============================================================================
// Metafunctions
// ============================================================================


// overrides for AligExtend general case
template <typename TSpec, typename TAlignConfig, typename TGapCosts,
          typename TTraceSwitch>
struct SetupAlignmentProfile_<AlignExtend_<TSpec>, TAlignConfig, TGapCosts,
                              TTraceSwitch>
{
    typedef DPProfile_<AlignExtend_<TSpec>, TGapCosts, TracebackOn<> > Type;
};

template <typename TSpec>
struct TraceTail_<AlignExtend_<TSpec> > : False
{};

template <typename TSpec>
struct TraceHead_<AlignExtend_<TSpec> > : True
{};

template <typename TSpec>
struct IsFreeEndGap_<AlignExtend_<TSpec>, DPLastRow> : True
{};

template <typename TSpec>
struct IsFreeEndGap_<AlignExtend_<TSpec>, DPLastColumn> : True
{};



// ----------------------------------------------------------------------------
// Class DPMetaColumn_                                             [FullColumn]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TGapCosts, typename TTraceback,
          typename TColumnType>
struct DPMetaColumn_<DPProfile_<AlignExtend_<TSpec>, TGapCosts,
                                TTraceback>,
                     MetaColumnDescriptor<TColumnType, FullColumn> >
{

    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionZero,
                        RecursionDirectionHorizontal
                       >::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionAll
                       >::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionAll
                       >::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};


// ----------------------------------------------------------------------------
// Class DPMetaColumn_                                       [PartialColumnTop]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TGapCosts, typename TTraceback,
          typename TColumnType>
struct DPMetaColumn_<DPProfile_<AlignExtend_<TSpec>, TGapCosts,
                                TTraceback>,
                     MetaColumnDescriptor<TColumnType, PartialColumnTop> >
{

    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionZero,
                        RecursionDirectionHorizontal
                       >::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionAll
                       >::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionLowerDiagonal
                       >::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn_                                    [PartialColumnMiddle]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TGapCosts, typename TTraceback,
          typename TColumnType>
struct DPMetaColumn_<DPProfile_<AlignExtend_<TSpec>, TGapCosts,
                                TTraceback>,
                     MetaColumnDescriptor<TColumnType, PartialColumnMiddle> >
{
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionZero,
                        RecursionDirectionUpperDiagonal
                       >::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionAll
                       >::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionLowerDiagonal
                       >::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn_                                    [PartialColumnBottom]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TGapCosts, typename TTraceback,
          typename TColumnType>
struct DPMetaColumn_<DPProfile_<AlignExtend_<TSpec>, TGapCosts,
                                TTraceback>,
                     MetaColumnDescriptor<TColumnType, PartialColumnBottom> >
{
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionZero,
                        RecursionDirectionUpperDiagonal
                       >::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionAll
                       >::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType,
                                   DPInitialColumn>, RecursionDirectionVertical,
                        RecursionDirectionAll
                       >::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

}
#endif
