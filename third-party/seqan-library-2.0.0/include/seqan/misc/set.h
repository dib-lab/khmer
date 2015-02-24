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

#ifndef SEQAN_HEADER_MISC_SET_H
#define SEQAN_HEADER_MISC_SET_H

#include <set>
#include <seqan/misc/base.h>

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    //
    //    VectorSet_
    //
    //////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////
    // forward declaration

    template <typename TElement, typename TSpec>
    class VectorSet_;


    //////////////////////////////////////////////////////////////////////////////
    // internal set meta-functions

    template <typename TElement>
    struct VectorSetKeySize_ {
        enum { VALUE = ValueSize< typename Key<TElement>::Type >::VALUE };
    };


    template <typename TSet>
    struct SetSetVector_ {
        typedef void* Type;
    };
    template <typename TElement, typename TSpec>
    struct SetSetVector_< VectorSet_<TElement, TSpec> > {
        typedef String<bool, TSpec> Type;
    };
    template <typename TSet>
    struct SetSetVector_<TSet const> {
        typedef typename SetSetVector_<TSet>::Type const Type;
    };



    template <typename TSet>
    struct SetObjVector_ {
        typedef void* Type;
    };
    template <typename TKey, typename TObject, typename TPairSpec, typename TSpec>
    struct SetObjVector_< VectorSet_<Pair<TKey, TObject, TPairSpec>, TSpec> > {
        typedef String<TObject, TSpec> Type;
    };
    template <typename TSet>
    struct SetObjVector_<TSet const> {
        typedef typename SetObjVector_<TSet>::Type const Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // VectorSet_ class

    template <
        typename TElement = char,
        typename TSpec = Alloc<> /*Array<VectorSetKeySize_<TKey>::VALUE>*/
    >
    class VectorSet_ {
    public:
        typedef typename SetSetVector_<VectorSet_>::Type        TSetVector;
        typedef typename SetObjVector_<VectorSet_>::Type        TObjVector;
        typedef typename Size<VectorSet_>::Type                TSize;

        TSetVector    vector;
        TObjVector    obj;
        TSize        size;

        VectorSet_():
            obj(), size(0)
        {
            _autoSize(*this);
            clear(*this);
        }

        VectorSet_(TSize _vectorSize):
            obj(), size(0)
        {
            resize(vector, _vectorSize);
            clear(*this);
        }

        template <typename TSet_>
        inline void _autoSize(TSet_ &) {}

        template <typename TElement_>
        inline void _autoSize(VectorSet_<TElement_, Alloc<> > &) {
            resize(vector, (unsigned)ValueSize<TElement_>::VALUE);
        }

        template <typename TKey_, typename TObject_, typename TSpec_>
        inline void _autoSize(VectorSet_<Pair<TKey_, TObject_, TSpec_>, Alloc<> > &) {
            resize(vector, (unsigned)ValueSize<TKey_>::VALUE);
            resize(obj, (unsigned)ValueSize<TKey_>::VALUE);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // set meta-functions

    template <typename TElement, typename TSpec>
    struct Value< VectorSet_<TElement, TSpec> > {
        typedef TElement Type;
    };
    template <typename TElement, typename TSpec>
    struct Size< VectorSet_<TElement, TSpec> >:
        Size< SetSetVector_< VectorSet_<TElement, TSpec> > > {};

    template <typename TElement, typename TSpec>
    struct Key< VectorSet_<TElement, TSpec> > :
        Key<TElement> {};

    template <typename TElement, typename TSpec>
    struct Object< VectorSet_<TElement, TSpec> > :
        Object<TElement> {};

    template <typename TObject, typename TSpec>
    inline typename Size< VectorSet_<TObject, TSpec> >::Type
    length(VectorSet_<TObject, TSpec> const &set) {
        return set.size;
    }



    //////////////////////////////////////////////////////////////////////////////
    //

    struct VectorSetIterator_;
    typedef Tag<VectorSetIterator_> VectorSetIterator__;

    template <typename TVectorSet>
    class Iter< TVectorSet, VectorSetIterator__ >
    {
        typedef typename SetSetVector_<TVectorSet>::Type        TSetVector;
        typedef typename SetObjVector_<TVectorSet>::Type        TObjVector;

        typedef Iter                                            iterator;
        typedef typename Value<TSetVector>::Type                TValue, value;
        typedef typename Iterator<TSetVector, Rooted>::Type        TSetIter;
        typedef typename Iterator<TObjVector, Standard>::Type    TObjIter;

    public:
        TSetIter    ptr;
        TObjIter    obj;

        Iter():
            obj(NULL) {}

        Iter(TSetIter _ptr, TObjIter _obj):
            ptr(_ptr), obj(_obj)
        {
            while (!atEnd(ptr) && !(*ptr)) {
                ++ptr;
                ++obj;
            }
        }

//        inline TValue operator*() { return TValue(ptr - begin, *obj); }

        inline iterator operator++() {
            ++ptr;
            ++obj;
            while (!atEnd(ptr) && !(*ptr)) {
                ++ptr;
                ++obj;
            }
            return *this;
        }
    };


    template <typename TObject, typename TSpec, typename TIteratorSpec>
    struct Iterator< VectorSet_<TObject, TSpec>, TIteratorSpec > {
        typedef Iter<VectorSet_<TObject, TSpec>, VectorSetIterator__> Type;
    };
    template <typename TObject, typename TSpec, typename TIteratorSpec>
    struct Iterator< VectorSet_<TObject, TSpec> const, TIteratorSpec> {
        typedef Iter<VectorSet_<TObject, TSpec> const, VectorSetIterator__> Type;
    };

    template <typename TValue, typename TSpec>
    finline void
    clear(VectorSet_<TValue, TSpec> &set) {
        arrayFill(begin(set.vector), end(set.vector), false);
        set.size = 0;
    }

    template <typename TElement, typename TSetKey, typename TSpec>
    finline void
    insert(TElement const &element, VectorSet_<TSetKey, TSpec> &set) {
        if (!set.vector[(unsigned)(element)]) {
            ++set.size;
            set.vector[(unsigned)(element)] = true;
        }
    }
    template <
        typename TElement,
        typename TSetKey,
        typename TSetObject,
        typename TPairSpec,
        typename TSpec>
    finline void
    insert(TElement const &element, VectorSet_< Pair<TSetKey, TSetObject, TPairSpec>, TSpec > &set) {
        if (!set.vector[(unsigned)(keyOf(element))]) {
            ++set.size;
            set.vector[(unsigned)(keyOf(element))] = true;
        }
        set.obj[(unsigned)(keyOf(element))] = objectOf(element);
    }

    template <typename TKey, typename TValue, typename TSpec>
    finline void
    erase(TKey const &key, VectorSet_<TValue, TSpec> &set) {
        if (set.vector[(unsigned)(key)]) {
            --set.size;
            set.vector[(unsigned)(key)] = false;
        }
    }

    template <typename TKey, typename TValue, typename TSpec>
    finline bool
    in(TKey const &key, VectorSet_<TValue, TSpec> const &set) {
        return set.vector[(unsigned)(key)];
    }

    template <typename TKey, typename TKey2, typename TSpec>
    inline typename Iterator< VectorSet_<TKey2, TSpec> >::Type
    find(TKey const &key, VectorSet_<TKey2, TSpec> &set) {
        if (in(key, set))
            return Iter<VectorSet_<TKey2, TSpec>, VectorSetIterator__>
                (begin(set.vector, Rooted()) + (unsigned)(key),
                 begin(set.obj, Standard()) + (unsigned)(key));
        else
            return end(set);
    }
    template <typename TKey, typename TKey2, typename TSpec>
    inline typename Iterator< VectorSet_<TKey2, TSpec> const>::Type
    find(TKey const &key, VectorSet_<TKey2, TSpec> const &set) {
        if (in(key, set))
            return Iter<VectorSet_<TKey2, TSpec> const, VectorSetIterator__>
                (begin(set.vector, Rooted()) + (unsigned)(key),
                 begin(set.obj, Standard()) + (unsigned)(key));
        else
            return end(set);
    }

    template <typename TElement, typename TSpec>
    inline typename Iterator< VectorSet_<TElement, TSpec> >::Type
    begin(VectorSet_<TElement, TSpec> &set) {
        return Iter<VectorSet_<TElement, TSpec>, VectorSetIterator__>
            (begin(set.vector, Rooted()), begin(set.obj, Standard()));
    }
    template <typename TElement, typename TSpec>
    inline typename Iterator< VectorSet_<TElement, TSpec> const>::Type
    begin(VectorSet_<TElement, TSpec> const &set) {
        return Iter<VectorSet_<TElement, TSpec> const, VectorSetIterator__>
            (begin(set.vector, Rooted()), begin(set.obj, Standard()));
    }

    template <typename TElement, typename TSpec>
    inline typename Iterator< VectorSet_<TElement, TSpec> >::Type
    end(VectorSet_<TElement, TSpec> &set) {
        return Iter<VectorSet_<TElement, TSpec>, VectorSetIterator__>
            (end(set.vector, Rooted()), begin(set.obj, Standard()));
    }
    template <typename TElement, typename TSpec>
    inline typename Iterator< VectorSet_<TElement, TSpec> const>::Type
    end(VectorSet_<TElement, TSpec> const &set) {
        return Iter<VectorSet_<TElement, TSpec> const, VectorSetIterator__>
            (end(set.vector, Rooted()), begin(set.obj, Standard()));
    }

    template <typename TSet>
    inline bool
    operator==(Iter<TSet, VectorSetIterator__> const &a, Iter<TSet, VectorSetIterator__> const &b) {
        return a.ptr == b.ptr;
    }
    template <typename TSet>
    inline bool
    operator!=(Iter<TSet, VectorSetIterator__> const &a, Iter<TSet, VectorSetIterator__> const &b) {
        return a.ptr != b.ptr;
    }
    template <typename TSet>
    inline bool
    eof(Iter<TSet, VectorSetIterator__> const &a) {
        return atEnd(a.ptr);
    }
    template <typename TSet>
    inline bool
    atEnd(Iter<TSet, VectorSetIterator__> const &a) {
        return atEnd(a.ptr);
    }

    //////////////////////////////////////////////////////////////////////////////

    template <typename TSet>
    inline typename Key<TSet>::Type
    keyOf(Iter<TSet, VectorSetIterator__> const &it) {
        return position(it.ptr);
    }

    template <typename TSet>
    inline typename Object<TSet>::Type &
    objectOf(Iter<TSet, VectorSetIterator__> const &it) {
        return *it.obj;
    }
    template <typename TSet>
    inline typename Object<TSet>::Type const &
    objectOf(Iter<TSet const, VectorSetIterator__> const &it) {
        return *it.obj;
    }



//____________________________________________________________________________




    //////////////////////////////////////////////////////////////////////////////
    //
    //    STL set adaptation
    //
    //////////////////////////////////////////////////////////////////////////////


    template <typename TElement>
    struct SetLess_ : public std::binary_function<TElement, TElement, bool>
    {
        // key less
        inline bool operator() (TElement const &a, TElement const &b) {
            return keyOf(a) < keyOf(b);
        }
    };

/* DISABLED: this part of code interferes with map_adapter_stl.h

    template <typename TElement>
    struct Value< std::set<TElement> > {
        typedef TElement Type;
    };
    template <typename TElement>
    struct Size< std::set<TElement> > {
        typedef typename std::set<TElement>::size_type Type;
    };

    template <typename TElement>
    struct Key< std::set<TElement> > :
        Key<TElement> {};

    template <typename TElement>
    struct Object< std::set<TElement> > :
        Object<TElement> {};



    template <typename TKey>
    inline void
    clear(std::set<TKey> &set) {
        set.clear();
    }

    template <typename TElement, typename TSetKey>
    inline void
    insert(TElement const &element, std::set<TSetKey> &set) {
        set.insert(element);
    }

    template <typename TKey, typename TSetKey>
    inline void
    erase(TKey const &key, std::set<TSetKey> &set) {
        set.erase(key);
    }
    template <typename TKey, typename TSetKey, typename TSetObject, typename TPairSpec>
    inline void
    erase(TKey const &key, std::set< Pair<TSetKey, TSetObject, TPairSpec> > &set) {
        set.erase(Pair<TSetKey, TSetObject, TPairSpec>(key, TSetObject()));
    }

    template <typename TKey, typename TSetKey>
    inline bool
    in(TKey const &key, std::set<TSetKey> const &set) {
        return set.count(key) != 0;
    }
    template <typename TKey, typename TSetKey, typename TSetObject, typename TPairSpec>
    inline bool
    in(TKey const &key, std::set<Pair<TSetKey, TSetObject, TPairSpec> > const &set) {
        return set.count(Pair<TSetKey, TSetObject, TPairSpec>(key, TSetObject())) != 0;
    }

    template <typename TKey>
    inline typename Size< std::set<TKey> >::Type
    length(std::set<TKey> const &set) {
        return set.size();
    }

    template <typename TKey, typename TSetKey>
    inline typename Iterator< std::set<TSetKey> >::Type
    find(TKey const &key, std::set<TSetKey> &set) {
        return set.find(key);
    }
    template <typename TKey, typename TSetKey>
    inline typename Iterator< std::set<TSetKey> const>::Type
    find(TKey const &key, std::set<TSetKey> const &set) {
        return set.find(key);
    }
    template <typename TKey, typename TSetKey, typename TSetObject, typename TPairSpec>
    inline typename Iterator< std::set<Pair<TSetKey, TSetObject, TPairSpec> > >::Type
    find(TKey const &key, std::set<Pair<TSetKey, TSetObject> > &set) {
        return set.find(Pair<TSetKey, TSetObject>(key, TSetObject()));
    }
    template <typename TKey, typename TSetKey, typename TSetObject>
    inline typename Iterator< std::set<Pair<TSetKey, TSetObject> > const>::Type
    find(TKey const &key, std::set<Pair<TSetKey, TSetObject> > const &set) {
        return set.find(Pair<TSetKey, TSetObject>(key, TSetObject()));
    }



    //////////////////////////////////////////////////////////////////////////////

    template <typename TObject>
    struct Iterator< std::set<TObject> > {
        typedef typename std::set<TObject>::iterator Type;
    };
    template <typename TObject>
    struct Iterator< std::set<TObject> const > {
        typedef typename std::set<TObject>::const_iterator Type;
    };


    template <typename TObject>
    typename Iterator< std::set<TObject> >::Type begin(std::set<TObject> &set) {
        return set.begin();
    }
    template <typename TObject>
    typename Iterator< std::set<TObject> const >::Type begin(std::set<TObject> const &set) {
        return set.begin();
    }
    template <typename TObject>
    typename Iterator< std::set<TObject> >::Type end(std::set<TObject> &set) {
        return set.end();
    }
    template <typename TObject>
    typename Iterator< std::set<TObject> const >::Type end(std::set<TObject> const &set) {
        return set.end();
    }

    //////////////////////////////////////////////////////////////////////////////

    template <typename TElement>
    inline typename Key< std::set<TElement> >::Type &
    keyOf(typename std::set<TElement>::iterator &it) {
        return keyOf(*it);
    }
    template <typename TElement>
    inline typename Key< std::set<TElement> >::Type &
    keyOf(typename std::set<TElement>::iterator const &it) {
        return keyOf(*it);
    }
    template <typename TElement>
    inline typename Key< std::set<TElement> >::Type const &
    keyOf(typename std::set<TElement>::const_iterator &it) {
        return keyOf(*it);
    }
    template <typename TElement>
    inline typename Key< std::set<TElement> >::Type const &
    keyOf(typename std::set<TElement>::const_iterator const &it) {
        return keyOf(*it);
    }

    template <typename TElement>
    inline typename Object< std::set<TElement> >::Type &
    objectOf(typename std::set<TElement>::iterator &it) {
        return objectOf(*it);
    }
    template <typename TElement>
    inline typename Object< std::set<TElement> >::Type &
    objectOf(typename std::set<TElement>::iterator const &it) {
        return objectOf(*it);
    }
    template <typename TElement>
    inline typename Object< std::set<TElement> >::Type const &
    objectOf(typename std::set<TElement>::const_iterator &it) {
        return objectOf(*it);
    }
    template <typename TElement>
    inline typename Object< std::set<TElement> >::Type const &
    objectOf(typename std::set<TElement>::const_iterator const &it) {
        return objectOf(*it);
    }


*/

//____________________________________________________________________________




    //////////////////////////////////////////////////////////////////////////////
    //
    //    Set chooser
    //
    //////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////
    // Set meta-function to choose an efficient implementation

    template <typename TKey>
    struct Set {
        typedef std::set<TKey> Type;
    };
    template <typename TKey, typename TObject, typename TPairSpec>
    struct Set< Pair<TKey, TObject, TPairSpec> > {
        typedef std::set<
            Pair<TKey, TObject, TPairSpec>,
            SetLess_< Pair<TKey, TObject, TPairSpec> > > Type;
    };


    template <typename TValue, typename TSpec>
    struct Set< SimpleType<TValue, TSpec> > {
        typedef VectorSet_< SimpleType<TValue, TSpec> > Type;
    };
    template <typename TValue, typename TSpec, typename TObject, typename TPairSpec>
    struct Set< Pair<SimpleType<TValue, TSpec>, TObject, TPairSpec> > {
        typedef VectorSet_< Pair<SimpleType<TValue, TSpec>, TObject, TPairSpec> > Type;
    };


    template <>
    struct Set<char> {
        typedef VectorSet_<char> Type;
    };
    template <typename TObject, typename TPairSpec>
    struct Set< Pair<char, TObject, TPairSpec> > {
        typedef VectorSet_< Pair<char, TObject, TPairSpec> > Type;
    };


    template <>
    struct Set<signed char> {
        typedef VectorSet_<signed char> Type;
    };
    template <typename TObject, typename TPairSpec>
    struct Set< Pair<signed char, TObject, TPairSpec> > {
        typedef VectorSet_< Pair<signed char, TObject, TPairSpec> > Type;
    };

    template <>
    struct Set<unsigned char> {
        typedef VectorSet_<unsigned char> Type;
    };
    template <typename TObject, typename TPairSpec>
    struct Set< Pair<unsigned char, TObject, TPairSpec> > {
        typedef VectorSet_< Pair<unsigned char, TObject, TPairSpec> > Type;
    };

}

#endif
