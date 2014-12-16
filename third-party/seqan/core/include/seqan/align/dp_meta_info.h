// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Renï¿½ Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Here are defined the strategies for the different alignment algorithms.
// All classes are only used to determine the correct computational state
// of a particular alignment algorithm depending on its profile at a
// particular time.
// All classes are only used on a meta-level.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_META_INFO_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_META_INFO_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DPInitialColumn
// ----------------------------------------------------------------------------

// Specifies the first column of the dp matrix.
struct DPInitialColumn_;
typedef Tag<DPInitialColumn_> DPInitialColumn;

// ----------------------------------------------------------------------------
// Tag DPInnerColumn
// ----------------------------------------------------------------------------

// Specifies any inner column of the dp matrix between the first and the last
// column.
struct DPInnerColumn_;
typedef Tag<DPInnerColumn_> DPInnerColumn;

// ----------------------------------------------------------------------------
// Tag DPFinalColumn
// ----------------------------------------------------------------------------

// Specifies the last column of a dp matrix.
struct DPFinalColumn_;
typedef Tag<DPFinalColumn_> DPFinalColumn;


// The TColumnProperty determines the property of the column (if it is initial, inner, or last column)
// The TLocation determines how the column is organized in the matrix.
// It can have the values: FullColumn, PartialColumnTop, PartialColumnMiddle, PartialColumnBottom.
template <typename TColumnProperty_, typename TLocation_>
struct MetaColumnDescriptor
{
    typedef TColumnProperty_ TColumnProperty;
    typedef TLocation_ TLocation;
};

// ----------------------------------------------------------------------------
// Tag FullColumn
// ----------------------------------------------------------------------------

// Columns that span over the complete dp-matrix. (Unbanded alignemnts)
struct FullColumn_;
typedef Tag<FullColumn_> FullColumn;

// ----------------------------------------------------------------------------
// Tag PartialColumnTop
// ----------------------------------------------------------------------------

// Columns that begin in the first row, but does not at the last row of the dp-matrix
struct PartialColumnTop_;
typedef Tag<PartialColumnTop_> PartialColumnTop;

// ----------------------------------------------------------------------------
// Tag PartialColumnMiddle
// ----------------------------------------------------------------------------

// Columns that are not attached to the begin or end of the dp-matrix.
struct PartialColumnMiddle_;
typedef Tag<PartialColumnMiddle_> PartialColumnMiddle;

// ----------------------------------------------------------------------------
// Tag PartialColumnBottom
// ----------------------------------------------------------------------------

// Columns that end in the last row of the dp-matrix but do not start in the first row.
struct PartialColumnBottom_;
typedef Tag<PartialColumnBottom_> PartialColumnBottom;


// The cell specifiers are used to determine the cell of a column.
// There are three different cell specifiers used for the first cell, the
// inner cell and the last cell of a column.

// ----------------------------------------------------------------------------
// Tag FirstCell
// ----------------------------------------------------------------------------

struct FirstCell_;
typedef Tag<FirstCell_> FirstCell;

// ----------------------------------------------------------------------------
// Tag InnerCell
// ----------------------------------------------------------------------------

struct InnerCell_;
typedef Tag<InnerCell_> InnerCell;

// ----------------------------------------------------------------------------
// Tag LastCell
// ----------------------------------------------------------------------------

struct LastCell_;
typedef Tag<LastCell_> LastCell;


// ----------------------------------------------------------------------------
// Class DPMetaCell_
// ----------------------------------------------------------------------------

// Keeps meta-information of a cell in the dp-matrix. It stores the recursion direction
// for a particular cell type and whether it can be tracked or not.

template <typename TDPRecursionDirection, typename TTrackFlag>
struct DPMetaCell_ {};

// ----------------------------------------------------------------------------
// Class DPMetaColumn_
// ----------------------------------------------------------------------------

// Keeps meta-information of an entire column in the dp-matrix. Depending on the chosen
// DPProfile and the current column descriptor it selects for the different column locations and
// cell specifiers the correct meta-information about how is which cell computed and which cell is
// tracked.

template <typename TDPProfile, typename TColumnDescriptor>
struct DPMetaColumn_ {};


// ----------------------------------------------------------------------------
// Class DPMetaColumn_											   [FullColumn]
// ----------------------------------------------------------------------------

template <typename TDPProfile, typename TColumnType>
struct DPMetaColumn_<TDPProfile, MetaColumnDescriptor<TColumnType, FullColumn> >
{
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> Horizontal | Zero, All, All
    // If FinalColumn -> Horizontal | Zero, All, All

    typedef typename If<Or<IsSameType<TColumnType, DPInitialColumn>,
                           IsFreeEndGap_<TDPProfile, DPFirstRow> >, RecursionDirectionZero, RecursionDirectionHorizontal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeLastCell_;

    // If Local
    // If InitialColumn -> True, True, True
    // If InnerColumn -> True, True, True
    // If FinalColumn -> True, True, True

    // If Global
    // If InitialColumn -> False, False, False | True (if DPLastRow)
    // If InnerColumn -> False, False, False | True (if DPLastRow)
    // If FinalColumn -> False | True, False | True (if DPLastColumn), True

    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >             // check this if he is really entitled to find the maximum here.
                           >, True, False>::Type TrackFlagFirstCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagInnerCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           Or<IsSameType<TColumnType, DPFinalColumn>,
                              IsFreeEndGap_<TDPProfile, DPLastRow> > >, True, False>::Type TrackFlagLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, TrackFlagFirstCell_> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, TrackFlagInnerCell_> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, TrackFlagLastCell_> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn_									     [PartialColumnTop]
// ----------------------------------------------------------------------------

template <typename TDPProfile, typename TColumnType>
struct DPMetaColumn_<TDPProfile, MetaColumnDescriptor<TColumnType, PartialColumnTop> >
{
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // How does the recursion directions look like?

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> Horizontal | Zero, All, LowerBand
    // If FinalColumn -> Horizontal | Zero, All, LowerBand

    typedef typename If<Or<IsSameType<TColumnType, DPInitialColumn>,
                           IsFreeEndGap_<TDPProfile, DPFirstRow> >, RecursionDirectionZero, RecursionDirectionHorizontal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionLowerDiagonal>::Type TRecursionTypeLastCell_;

    // If Local
    // If InitialColumn -> True, True, True
    // If InnerColumn -> True, True, True
    // If FinalColumn -> True, True, True

    // If Global
    // If InitialColumn -> False, False, False
    // If InnerColumn -> False, False, False
    // If FinalColumn -> False | True, False | True, False | True (if DPLastColumn True)

    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagFirstCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagInnerCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, TrackFlagFirstCell_> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, TrackFlagInnerCell_> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, TrackFlagLastCell_> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn_									  [PartialColumnMiddle]
// ----------------------------------------------------------------------------

template <typename TDPProfile, typename TColumnType>
struct DPMetaColumn_<TDPProfile, MetaColumnDescriptor<TColumnType, PartialColumnMiddle> >
{
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> UpperDiagonal, All, LowerDiagonal
    // If FinalColumn -> UpperDiagonal, All, LowerDiagonal

    typedef typename If<IsSameType<TColumnType, DPInitialColumn>, RecursionDirectionZero, RecursionDirectionUpperDiagonal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionLowerDiagonal>::Type TRecursionTypeLastCell_;

    // If Local
    // If InitialColumn -> True, True, True
    // If InnerColumn -> True, True, True
    // If FinalColumn -> True, True, True

    // If Global
    // If InitialColumn -> False, False, False
    // If InnerColumn -> False, False, False
    // If FinalColumn -> False | True, False | True, False | True (if DPLastColumn True)

    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagFirstCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagInnerCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, TrackFlagFirstCell_> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, TrackFlagInnerCell_> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, TrackFlagLastCell_> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn_									  [PartialColumnBottom]
// ----------------------------------------------------------------------------

template <typename TDPProfile, typename TColumnType>
struct DPMetaColumn_<TDPProfile, MetaColumnDescriptor<TColumnType, PartialColumnBottom> >
{
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> UpperDiagonal, All, All
    // If FinalColumn -> UpperDiagonal, All, All

    typedef typename If<IsSameType<TColumnType, DPInitialColumn>, RecursionDirectionZero, RecursionDirectionUpperDiagonal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeLastCell_;

    // If Local
    // If InitialColumn -> True, True, True
    // If InnerColumn -> True, True, True
    // If FinalColumn -> True, True, True

    // If Global
    // If InitialColumn -> False, False, False | True (if DPLastRow)
    // If InnerColumn -> False, False, False | True (if DPLastRow)
    // If FinalColumn -> False | True, False | True (if DPLastColumn), True (last is always true)

    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagFirstCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           And<IsSameType<TColumnType, DPFinalColumn>,
                               IsFreeEndGap_<TDPProfile, DPLastColumn> >
                           >, True, False>::Type TrackFlagInnerCell_;
    typedef typename If<Or<IsLocalAlignment_<TDPProfile>,
                           Or<IsSameType<TColumnType, DPFinalColumn>,
                              IsFreeEndGap_<TDPProfile, DPLastRow> > >, True, False>::Type TrackFlagLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, TrackFlagFirstCell_> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, TrackFlagInnerCell_> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, TrackFlagLastCell_> TLastCell_;
};


// ----------------------------------------------------------------------------
// Metafunction GetRecursionDirection_
// ----------------------------------------------------------------------------

// Returns the type of recursion for a given DPMetaCell object.
template <typename TMetaCell>
struct GetRecursionDirection_
{
    typedef Nothing Type;
};

template <typename TRecursionDirection, typename TTraceFlag>
struct GetRecursionDirection_<DPMetaCell_<TRecursionDirection, TTraceFlag> >
{
    typedef TRecursionDirection Type;
};

// ----------------------------------------------------------------------------
// Metafunction RecursionDirection_
// ----------------------------------------------------------------------------

// Returns the type of recursion for a given DPMetaColumn object and a given cell specifier.
template <typename TDPMetaColumn, typename TCellDescriptor>
struct RecursionDirection_ {};

template <typename TDPMetaColumn>
struct RecursionDirection_<TDPMetaColumn, FirstCell>
{
    typedef typename GetRecursionDirection_<typename TDPMetaColumn::TFirstCell_>::Type Type;
};

template <typename TDPMetaColumn>
struct RecursionDirection_<TDPMetaColumn, InnerCell>
{
    typedef typename GetRecursionDirection_<typename TDPMetaColumn::TInnerCell_>::Type Type;
};

template <typename TDPMetaColumn>
struct RecursionDirection_<TDPMetaColumn, LastCell>
{
    typedef typename GetRecursionDirection_<typename TDPMetaColumn::TLastCell_>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsTrackingEnabled_
// ----------------------------------------------------------------------------

// Returns an object that evaluates to true if for a given cell description the
// tracking was enabled. Otherwise the object evaluates to false.
template <typename TMetaCell>
struct IsTrackingEnabled_ :
    False {};

template <typename TRecursionDirection>
struct IsTrackingEnabled_<DPMetaCell_<TRecursionDirection, True> >:
    True {};


template <typename TDPMetaColumn, typename TCellDescriptor>
struct TrackingEnabled_ :
    False {};

template <typename TDPMetaColumn>
struct TrackingEnabled_<TDPMetaColumn, FirstCell>:
    IsTrackingEnabled_<typename TDPMetaColumn::TFirstCell_>{};

template <typename TDPMetaColumn>
struct TrackingEnabled_<TDPMetaColumn, InnerCell>:
    IsTrackingEnabled_<typename TDPMetaColumn::TInnerCell_>{};

template <typename TDPMetaColumn>
struct TrackingEnabled_<TDPMetaColumn, LastCell>:
    IsTrackingEnabled_<typename TDPMetaColumn::TLastCell_>{};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_META_INFO_H_
