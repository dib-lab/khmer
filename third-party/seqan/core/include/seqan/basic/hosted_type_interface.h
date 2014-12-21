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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Interface specification functions for types that have a host.
// ==========================================================================

// TODO(holtgrew): We could add a HostedTypeConcept and make this a submodule of basic, e.g. basic/hosted.
// TODO(holtgrew): This looks a bit unused/underused.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Concept.HostedConcept Type
..summary:Concept for types that have a host.
..remarks:The functions of this concept assume that the hosted object exports a function $_dataHost$ that returns a reference to a holder type of $Host<T>::Type &$.

.Metafunction.Host.concept:Concept.HostedConcept Type
*/

// ============================================================================
// Metafunctions
// ============================================================================

/**
.Metafunction.Host
..cat:Basic
..summary:Type of the object a given object depends on.
..signature:Host<T>::Type
..param.T:Type for which the host type is determined.
..returns.param.Type:Host type of $T$.
..include:seqan/basic.h
*/

template <typename T>
struct Host;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function emptyHost()
// ----------------------------------------------------------------------------

/**
.Function.emptyHost
..cat:Dependent Object
..summary:Query emptiness state of a hosted object.
..signature:emptyHost(object)
..param.object:The object query state of host of.
..returns:$bool$, $true$ if the host is empty, $false$ otherwise.
..see:Function.empty
..concept:Concept.HostedConcept Type
..include:seqan/basic.h
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

/**
.Function.dependentHost
..cat:Dependent Object
..summary:Query dependent state of a hosted object.
..signature:clearHost(object)
..param.object:The object query state of host of.
..returns:$bool$, $true$ if the host is dependent, $false$ otherwise.
..see:Function.dependent
..concept:Concept.HostedConcept Type
..include:seqan/basic.h
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

/**
.Function.clearHost
..cat:Dependent Object
..summary:Clear the host of the given object.
..signature:clearHost(object)
..param.object:The object to clear the host of.
..see:Function.clear
..concept:Concept.HostedConcept Type
..include:seqan/basic.h
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

/**
.Function.createHost
..cat:Dependent Object
..summary:Construct the host of the given object.
..signature:createHost(object[, host])
..param.object:The object to copy construct the host of.
..param.host:The object to copy in host creation.
...type:nolink:$Host<T>::Type const &$
..remarks:If $host$ is given then it is used for copy creation.  Otherwise, the default constructor is used.
..see:Function.create
..concept:Concept.HostedConcept Type
..include:seqan/basic.h
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

/// TODO(holtgrew): Move documentation here?

///.Function.host.concept:Concept.HostedConcept Type

template <typename T>
inline typename Host<T>::Type &
host(T & me)
{
    SEQAN_CHECKPOINT;
    return value(_dataHost(me));
}

// TODO(holtgrew): Is this function unnecessary? Should be since the above one is catch-all.
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

/// TODO(holtgrew): Move documentation here?

///.Function.setHost.param.object.type:nolink:$Host<T>::Type &$
///.Function.setHost.concept:Concept.HostedConcept Type

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

/**
.Function.assignHost
..cat:Dependent Object
..summary:Assign to the host of a given value.
..signature:assignHost(object, host)
..param.object:The object to assign the host of.
..param.host:The object to assign as host.
...type:nolink:$Host<T>::Type const &$
..see:Function.assign
..concept:Concept.HostedConcept Type
..include:seqan/basic.h
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

/**
.Function.moveHost
..cat:Dependent Object
..summary:Assign to the host of a given value.
..signature:assignHost(object, host)
..param.object:The object to move-assign the host of.
..param.host:The object to move-assign as host.
...type:nolink:$Host<T>::Type &$
..see:Function.move
..concept:Concept.HostedConcept Type
..include:seqan/basic.h
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

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_
