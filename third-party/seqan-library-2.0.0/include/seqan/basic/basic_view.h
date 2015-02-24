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

#ifndef SEQAN_BASIC_VIEW_H
#define SEQAN_BASIC_VIEW_H

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

/*!
 * @mfn View
 * @headerfile <seqan/basic.h>
 * @brief Converts a given type into its view type.
 *
 * @signature View<TObject>::Type;
 *
 * @tparam TObject The type to be converted into a view type.
 * @return Type The resulting view type.
 *
 * This metafunction is used to convert device containers into views of device containers.
 * Subsequently, the view of a device container can be safely passed to and used in device space.
 * On the host, a view of a @link String @endlink is equivalent to an @link SegmentableConcept#Infix @endlink of the complete string.
 * @link RemoveView @endlink is the inverse of this metafunction.
 *
 * @see RemoveView
 * @see Device
 * @see SegmentableConcept#Infix
 */

template <typename TObject>
struct View
{
    typedef TObject Type;
};

template <typename TObject>
struct View<TObject const>
{
    typedef typename View<TObject>::Type const  Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView
// ----------------------------------------------------------------------------

/*!
 * @mfn RemoveView
 * @headerfile <seqan/basic.h>
 * @brief Converts a given view type into its original type.
 *
 * @signature RemoveView<TObject>::Type;
 *
 * @tparam TObject The view type to be converted into its original type.
 * @return Type The resulting original type.
 *
 * @link View @endlink is the inverse of this metafunction.
 *
 * @see View
 */

template <typename TObject>
struct RemoveView
{
    typedef TObject Type;
};

template <typename TObject>
struct RemoveView<TObject const>
{
    typedef typename RemoveView<TObject>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView
// ----------------------------------------------------------------------------

/*!
 * @mfn IsView
 * @headerfile <seqan/basic.h>
 * @brief Tests if a given type is a view type.
 *
 * @signature IsView<TObject>::Type;
 *
 * @tparam TObject The type to be tested for being a view type.
 * @return Type @link LogicalValuesTags#True @endlink or @link LogicalValuesTags#False @endlink.
 *
 * @see View
 * @see RemoveView
 */

template <typename TObject>
struct IsView : public False {};

template <typename TObject>
struct IsView<TObject const> : public IsView<TObject> {};

// ----------------------------------------------------------------------------
// Metafunction IfView
// ----------------------------------------------------------------------------

template <typename TObject, typename T1, typename T2>
struct IfView
{
    typedef typename If<IsView<TObject>, T1, T2>::Type  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

/*!
 * @fn TView view
 * @headerfile <seqan/basic.h>
 * @brief Returns the view of a given object.
 *
 * @signature TView view(object);
 *
 * @param[in] object A generic object.
 * @return TView The @link View @endlink type of the given object.
 *
 * @see View
 * @see IsView
 */

template <typename TObject>
inline typename View<TObject>::Type
view(TObject & object)
{
    return typename View<TObject>::Type(object);
}

template <typename TObject>
inline typename View<TObject const>::Type
view(TObject const & object)
{
    return typename View<TObject const>::Type(object);
}

template <typename TObject>
inline typename View<TObject>::Type
view(TObject * object)
{
    return typename View<TObject>::Type(value(object));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_VIEW_H
