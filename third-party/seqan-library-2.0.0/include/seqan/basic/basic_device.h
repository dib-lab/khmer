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

#ifndef SEQAN_BASIC_DEVICE_H
#define SEQAN_BASIC_DEVICE_H

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Execution space tags
// ----------------------------------------------------------------------------

struct ExecHost_;
struct ExecDevice_;

typedef Tag<ExecHost_>   ExecHost;
typedef Tag<ExecDevice_> ExecDevice;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

/*!
 * @mfn Device
 * @headerfile <seqan/basic.h>
 * @brief Converts a given type into one that lives on a device.
 *
 * @signature Device<TObject>::Type;
 * @tparam TObject The type to be converted into a device type.
 * @return Type The resulting device type.
 *
 * This metafunction is used to convert host containers into device containers.
 *
 * @see View
 */

template <typename TObject = void>
struct Device
{
    typedef TObject Type;
};

template <typename TObject>
struct Device<TObject const>
{
    typedef typename Device<TObject>::Type const    Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsDevice
// ----------------------------------------------------------------------------

/*!
 * @mfn IsDevice
 * @headerfile <seqan/basic.h>
 * @brief Tests if a given type is a device type.
 *
 * @signature IsDevice<TObject>::Type;
 * @tparam TObject The type to be tested for being a device type.
 * @return Type @link LogicalValuesTags#True @endlink or @link LogicalValuesTags#False @endlink.
 *
 * @see Device
 */

template <typename TObject>
struct IsDevice : public False {};

template <typename TObject>
struct IsDevice<TObject const> : public IsDevice<TObject> {};

// ----------------------------------------------------------------------------
// Metafunction IfDevice
// ----------------------------------------------------------------------------

template <typename TObject, typename T1, typename T2>
struct IfDevice
{
    typedef typename If<IsDevice<TObject>, T1, T2>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction ExecSpace
// ----------------------------------------------------------------------------

template <typename TObject>
struct ExecSpace
{
    typedef typename If<IsDevice<TObject>, ExecDevice, ExecHost>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction ExecSpec
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec = void>
struct ExecSpec
{
    typedef typename IfDevice<TObject, Device<TSpec>, TSpec>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction CtaSize
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec = void>
struct CtaSize
{
    static const unsigned VALUE = 256;
};

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_DEVICE_H
