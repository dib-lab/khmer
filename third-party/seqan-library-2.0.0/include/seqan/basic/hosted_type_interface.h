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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Interface specification functions for types that have a host.
// ==========================================================================

// TODO(holtgrew): We could add a HostedTypeConcept and make this a submodule of basic, e.g. basic/hosted.
// TODO(holtgrew): This looks a bit unused/underused.

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @concept HostedConcept
 * @brief Concept for types that have a host.
 *
 * @signature concept HostedConcept;
 *
 * @section Remarks
 *
 * The functions of this concept assume that the hosted object exports a function <tt>_dataHost</tt> that returns a
 * reference to a holder type of <tt>Host&lt;T&gt;::Type &amp;</tt>.
 */

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn HostedConcept#Host
 * @brief Type of the object a given object depends on.
 *
 * @signature Host<T>::Type
 *
 * @tparam T Type for which the host type is determined.
 * @return Type The Host type.
 */

template <typename T>
struct Host;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function emptyHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#emptyHost
 * @brief Query emptiness state of a hosted object.
 *
 * @signature bool emptyHost(object);
 *
 * @param[in] object The object query state of host of.
 *
 * @return bool <tt>true</tt> if the host is empty, <tt>false</tt> otherwise.
 */

template <typename T>
inline bool
emptyHost(T const & me)
{
    SEQAN_CHECKPOINT;
    return empty(_dataHost(me));
}

// ----------------------------------------------------------------------------
// Function dependentHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#dependentHost
 * @brief Query dependent state of a hosted object.
 *
 * @signature void clearHost(object);
 *
 * @param[in] object The object query state of host of.
 *
 * @return bool <tt>true</tt> if the host is dependent, <tt>false</tt> otherwise.
 */

template <typename T>
inline bool
dependentHost(T const & me)
{
    SEQAN_CHECKPOINT;
    return dependent(_dataHost(me));
}

// ----------------------------------------------------------------------------
// Function clearHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#clearHost
 * @brief Clear the host of the given object.
 *
 * @signature void clearHost(object);
 *
 * @param[in,out] object The object to clear the host of.
 */

template <typename T>
inline void
clearHost(T & me)
{
    SEQAN_CHECKPOINT;
    clear(_dataHost(me));
}

// ----------------------------------------------------------------------------
// Function createHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#createHost
 * @brief Construct the host of the given object.
 *
 * @signature void createHost(object[, host]);
 *
 * @param[in,out] object The object to copy construct the host of.
 * @param[in]     host   The object to copy in host creation.
 *
 * @section Remarks
 *
 * If <tt>host</tt> is given then it is used for copy creation.  Otherwise, the default constructor is used.
 */

template <typename T>
inline void
createHost(T & me)
{
    SEQAN_CHECKPOINT;
    create(_dataHost(me));
}

template <typename T, typename THost>
inline void
createHost(T & me,
           THost const & host_)
{
    SEQAN_CHECKPOINT;
    create(_dataHost(me), host_);
}

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#host
 * @brief The object a given object depends on.
 *
 * @signature THostRef host(object);
 *
 * @param[in] object An object.
 *
 * @return THostRef Reference to the host object.
 */

/// TODO(holtgrew): Move documentation here?

template <typename T>
inline typename Host<T>::Type &
host(T & me)
{
    SEQAN_CHECKPOINT;
    return value(_dataHost(me));
}

// TODO(holtgrew): Is this function unnecessary? Should be since the above one is catch-all.
// (weese:) No, the above wouldn't catch const refs.
template <typename T>
inline typename Host<T const>::Type &
host(T const & me)
{
    SEQAN_CHECKPOINT;
    return value(_dataHost(me));
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#setHost
 * @brief Sets the host of an object.
 *
 * @signature void setHost(object, host);
 *
 * @param[in,out] host   The new host. Types: String
 * @param[in]     object The object that will get a new host.
 *
 * @section Remarks
 *
 * After this operation, <tt>object</tt> depends on <tt>host</tt>.
 *
 * Note that setting the host can invalidate <tt>object</tt>.  For example, if one changes the host of a Segment object,
 * it is possible that begin- and end-position of the segment does not fit into the new host sequence.
 */

/// TODO(holtgrew): Move documentation here?

template <typename T, typename THost>
inline void
setHost(T & me,
        THost & host_)
{
    SEQAN_CHECKPOINT;
    setValue(_dataHost(me), host_);
}

template <typename T, typename THost>
inline void
setHost(T & me,
        THost const & host_)
{
    SEQAN_CHECKPOINT;
    setValue(_dataHost(me), host_);
}

// ----------------------------------------------------------------------------
// Function assignHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#assignHost
 * @brief Assign to the host of a given value.
 *
 * @signature void assignHost(object, host);
 *
 * @param[in,out] host   The object to assign as host.
 * @param[in]     object The object to assign the host of.
 */

template <typename T, typename THost>
inline void
assignHost(T & me,
           THost const & host_)
{
    SEQAN_CHECKPOINT;
    assignValue(_dataHost(me), host_);
}

// ----------------------------------------------------------------------------
// Function moveHost()
// ----------------------------------------------------------------------------

/*!
 * @fn HostedConcept#moveHost
 * @brief Move to the host of a given value.
 *
 * @signature void moveHost(object, host);
 *
 * @param[in,out] host   The object to move-assign as host.
 * @param[in,out] object The object to move-assign the host of.
 */

template <typename T, typename THost>
inline void
moveHost(T & me,
         THost & host_)
{
    SEQAN_CHECKPOINT;
    moveValue(_dataHost(me), host_);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_
