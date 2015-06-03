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

#ifndef SEQAN_INDEX_ITER_VSTREE_FACTORY_H
#define SEQAN_INDEX_ITER_VSTREE_FACTORY_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Class Factory
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec = void>
struct Factory;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
struct View<Factory<TObject, TSpec> >
{
    typedef Factory<typename View<TObject>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
struct IsView<Factory<TObject, TSpec> > : IsView<TObject> {};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Host<Factory<Iter<TIndex, VSTree<TSpec> > > >
{
    typedef TIndex  Type;
};

template <typename TIndex, typename TSpec>
struct Host<Factory<Iter<TIndex, VSTree<TSpec> > > const>
{
    typedef TIndex const    Type;
};

// ----------------------------------------------------------------------------
// Member Index_
// ----------------------------------------------------------------------------

struct Index_;

template <typename TIndex, typename TSpec>
struct Member<Factory<Iter<TIndex, VSTree<TSpec> > >, Index_>
{
    typedef typename IfView<TIndex, TIndex, Holder<TIndex> >::Type  Type;
};

// ----------------------------------------------------------------------------
// Member History_
// ----------------------------------------------------------------------------

struct History_;

template <typename TIndex, typename TSpec>
struct Member<Factory<Iter<TIndex, TSpec> >, History_>
{
    typedef typename HistoryStack_<Iter<TIndex, TSpec> >::Type   Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Factory
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Factory<Iter<TIndex, VSTree<TSpec> > >
{
    typename Member<Factory, Index_>::Type  _index;

    Factory() {}

    Factory(TIndex & index) :
        _index(index)
    {}
};

// ----------------------------------------------------------------------------
// Class Factory
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > :
    Factory<Iter<TIndex, VSTree<TSpec> > >
{
    typedef Factory<Iter<TIndex, VSTree<TSpec> > >  TBase;

    typename Member<Factory, History_>::Type    _history;
    typename Size<TIndex>::Type                 _maxHistoryLength;
    typename Size<Factory>::Type                _maxObjects;

    Factory() :
        TBase(),
        _maxHistoryLength(0),
        _maxObjects(0)
    {}

    Factory(TIndex & index) :
        TBase(index)
    {}

    template <typename TSize, typename THistorySize>
    Factory(TIndex & index, TSize maxObjects, THistorySize maxHistoryLength) :
        TBase(index)
    {
        setMaxObjects(*this, maxObjects);
        setMaxHistoryLength(*this, maxHistoryLength);
        build(*this);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline typename Host<Factory<Iter<TIndex, VSTree<TSpec> > > >::Type &
host(Factory<Iter<TIndex, VSTree<TSpec> > > & factory)
{
    return _host(factory, typename IsView<TIndex>::Type());
}

template <typename TIndex, typename TSpec>
inline SEQAN_HOST_DEVICE typename Host<Factory<Iter<TIndex, VSTree<TSpec> > > >::Type &
_host(Factory<Iter<TIndex, VSTree<TSpec> > > & factory, True const & /* isView */)
{
    return factory._index;
}

template <typename TIndex, typename TSpec>
inline typename Host<Factory<Iter<TIndex, VSTree<TSpec> > > >::Type &
_host(Factory<Iter<TIndex, VSTree<TSpec> > > & factory, False const & /* isView */)
{
    return value(factory._index);
}

template <typename TIndex, typename TSpec>
inline typename Host<Factory<Iter<TIndex, VSTree<TSpec> > > const>::Type &
host(Factory<Iter<TIndex, VSTree<TSpec> > > const & factory)
{
    return _host(factory, typename IsView<TIndex>::Type());
}

template <typename TIndex, typename TSpec>
inline SEQAN_HOST_DEVICE typename Host<Factory<Iter<TIndex, VSTree<TSpec> > > const>::Type &
_host(Factory<Iter<TIndex, VSTree<TSpec> > > const & factory, True const & /* isView */)
{
    return factory._index;
}

template <typename TIndex, typename TSpec>
inline typename Host<Factory<Iter<TIndex, VSTree<TSpec> > > const>::Type &
_host(Factory<Iter<TIndex, VSTree<TSpec> > > const & factory, False const & /* isView */)
{
    return value(factory._index);
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline typename View<Factory<Iter<TIndex, VSTree<TSpec> > > >::Type
view(Factory<Iter<TIndex, VSTree<TSpec> > > & factory)
{
    typename View<Factory<Iter<TIndex, VSTree<TSpec> > > >::Type factoryView;

    host(factoryView) = view(_host(factory, typename IsView<TIndex>::Type()));

    return factoryView;
}

template <typename TIndex, typename TSpec>
inline typename View<Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > >::Type
view(Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > & factory)
{
    typename View<Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > >::Type factoryView;

    host(factoryView) = view(_host(factory, typename IsView<TIndex>::Type()));
    factoryView._history = view(factory._history);
    factoryView._maxHistoryLength = factory._maxHistoryLength;
    factoryView._maxObjects = factory._maxObjects;

    return factoryView;
}

// ----------------------------------------------------------------------------
// Function setMaxHistoryLength()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TSize>
inline void
setMaxHistoryLength(Factory<Iter<TIndex, VSTree<TSpec> > > & /* factory */, TSize /* maxHistoryLength */)
{}

template <typename TIndex, typename TSpec, typename TSize>
inline void
setMaxHistoryLength(Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > & factory, TSize maxHistoryLength)
{
    factory._maxHistoryLength = maxHistoryLength;
}

// ----------------------------------------------------------------------------
// Function setMaxObjects()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TSize>
inline void
setMaxObjects(Factory<Iter<TIndex, VSTree<TSpec> > > & /* factory */, TSize /* maxObjects */)
{}

template <typename TIndex, typename TSpec, typename TSize>
inline void
setMaxObjects(Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > & factory, TSize maxObjects)
{
    factory._maxObjects = maxObjects;
}

// ----------------------------------------------------------------------------
// Function build()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline void
build(Factory<Iter<TIndex, VSTree<TSpec> > > & /* factory */)
{}

template <typename TIndex, typename TSpec>
inline void
build(Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > & factory)
{
    resize(factory._history, factory._maxObjects * factory._maxHistoryLength, Exact());
}

// ----------------------------------------------------------------------------
// Function getObject()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TId>
inline SEQAN_HOST_DEVICE Iter<TIndex, VSTree<TSpec> >
getObject(Factory<Iter<TIndex, VSTree<TSpec> > > & factory, TId /* objectId */)
{
    return Iter<TIndex, VSTree<TSpec> >(_host(factory, typename IsView<TIndex>::Type()));
}

// ----------------------------------------------------------------------------
// Function getObject()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TId>
inline SEQAN_HOST_DEVICE Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > >
getObject(Factory<Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > > & factory, TId objectId)
{
    SEQAN_ASSERT_LT(objectId, factory._maxObjects);

    Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > it(_host(factory, typename IsView<TIndex>::Type()));

    it.history._begin = begin(factory._history, Standard()) + objectId * factory._maxHistoryLength;
    it.history._end = it.history._begin;
    it.history._capacity = factory._maxHistoryLength;

    return it;
}

}

#endif  // #ifndef SEQAN_INDEX_ITER_VSTREE_FACTORY_H
