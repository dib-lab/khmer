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

#ifndef SEQAN_HEADER_MISC_SKIPLIST_H
#define SEQAN_HEADER_MISC_SKIPLIST_H


#include <seqan/random.h>


namespace SEQAN_NAMESPACE_MAIN
{

/*!
 * @class Skiplist
 * @extends Map
 * @headerfile <seqan/map.h>
 *
 * @brief General purpose map container.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class Map<TValue, Skiplist<TSpec> >;
 *
 * @tparam TSpec  The specializing type.
 * @tparam TValue The type of value stored in the map.
 *
 * The skiplist takes in average an oberhead of only two pointers per value stored in the map.
 */

//////////////////////////////////////////////////////////////////////////////
// Skiplist
//////////////////////////////////////////////////////////////////////////////

//forwards

template <typename TValue, typename TSpec>
class Map;

template <typename TValue, typename TSpec>
class SkiplistElement;

template <typename TValue, typename TSpec>
class SkiplistNext;

template <typename TValue, typename TSpec>
class SkiplistPath;

//////////////////////////////////////////////////////////////////////////////
// Tags

struct SkiplistIterator;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct AllocatorType;

template <typename TValue, typename TSpec>
struct AllocatorType<Map<TValue, Skiplist<TSpec> > >
{
    typedef Allocator<SimpleAlloc<> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value<Map<TValue, Skiplist<TSpec> > >
{
    typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct SkiplistElement_
{
    typedef char Type; //dummy implementation for VC++ bug
};

template <typename TValue, typename TSpec>
struct SkiplistElement_<Map<TValue, Skiplist<TSpec> > >
{
    typedef SkiplistElement<TValue, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Key<Map<TValue, Skiplist<TSpec> > >:
    Key<typename Value< Map<TValue, Skiplist<TSpec> > >::Type>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Cargo<Map<TValue, Skiplist<TSpec> > >:
    Cargo<typename Value< Map<TValue, Skiplist<TSpec> > >::Type>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator<Map<TValue, Skiplist<TSpec> >, TIteratorSpec >
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist_;
    typedef Iter<TSkiplist_, SkiplistIterator> Type;
};


//////////////////////////////////////////////////////////////////////////////
// Skiplist Map
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
class Map<TValue, Skiplist<TSpec> >
{
public:
    typedef typename AllocatorType<Map>::Type TAllocator;
    typedef SkiplistElement<TValue, TSpec> TElement;
    typedef typename Size<Map>::Type TSize;

    enum
    {
        MAX_HEIGHT = 28, //heights are in {0, 1, ..., MAX_HEIGHT-1}
        BLOCK_SIZE_1_ = sizeof(TElement) * 20, //store min. 20 elements
        BLOCK_SIZE_2_ = 0x200,    //minimal block size
        BLOCK_SIZE = (BLOCK_SIZE_1_ > BLOCK_SIZE_2_) ? BLOCK_SIZE_1_ : BLOCK_SIZE_2_ //block size is the max out of both values
    };

    Holder<TAllocator> data_allocator;
    TElement * data_recycle[MAX_HEIGHT];
    unsigned char * data_mem_begin;
    unsigned char * data_mem_end;

    SkiplistElement<TValue, TSpec> data_border;
    TSize data_length;
    unsigned char data_height;

    Rng<> rng;

    Map()
        : data_mem_begin(0)
        , data_mem_end(0)
        , data_length(0)
        , data_height(0)
        , rng(0)
    {
        for (unsigned char i = 0; i < MAX_HEIGHT; ++i)
        {
            data_recycle[i] = 0;
            valueConstruct(data_border.data_next + i, NonMinimalCtor());
        }
    }
    Map(Map const & other)
        : data_mem_begin(0)
        , data_mem_end(0)
        , data_length(0)
        , data_height(0)
        , rng(0)
    {
        assign(*this, other);
    }

    Map &
    operator=(Map const & other)
    {
        assign(*this, other);
        return *this;
    }

    template <typename TKey>
    inline typename MapValue<Map>::Type
    operator [] (TKey const & key)
    {
        return mapValue(*this, key);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
class SkiplistElement
{
public:
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef SkiplistNext<TValue, TSpec> TNext;

    enum
    {
        MAX_HEIGHT = TSkiplist::MAX_HEIGHT
    };

    TValue data_value;
    TNext data_next[MAX_HEIGHT]; //note: only parts of this array available

    //indirect constructor in _skiplistConstructElement
    //indirect destructor in _skiplistDestructElement
};

//////////////////////////////////////////////////////////////////////////////

//representation of "horizontal" pointer in skiplist
//can be overloaded if label is needed
template <typename TValue, typename TSpec>
class SkiplistNext
{
public:
    typedef SkiplistElement<TValue, TSpec> TElement;

    TElement * data_element;

    SkiplistNext() : data_element(0)
    {}

    SkiplistNext(NonMinimalCtor) : data_element(0)
    {}

    SkiplistNext(SkiplistNext const & other) : data_element(other.data_element)
    {}
};

//////////////////////////////////////////////////////////////////////////////

//represents a path from the root to the predecessor of a skiplist element
template <typename TValue, typename TSpec>
class SkiplistPath
{
public:
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef SkiplistElement<TValue, TSpec> TElement;

    enum
    {
        MAX_HEIGHT = TSkiplist::MAX_HEIGHT
    };

    TElement * data_elements[MAX_HEIGHT];
};

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
inline void
assign(Map<TValue, Skiplist<TSpec> > & target,
       Map<TValue, Skiplist<TSpec> > const & source)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef SkiplistPath<TValue, TSpec> TPath;
    typedef SkiplistElement<TValue, TSpec> TElement;
    typedef typename Iterator<TSkiplist>::Type TIterator;

    clear(target);

    TPath path;
    path.data_elements[0] = & target.data_border;

    for (TIterator it(source); !atEnd(it); goNext(it))
    {
        unsigned char height = _skiplistCreateHeight(target, path);
        TElement & el = _skiplistConstructElement(target, height, value(it));

        for (int i = 0; i <= height; ++i)
        {
            el.data_next[i].data_element = 0;
            path.data_elements[i]->data_next[i].data_element = & el;
            path.data_elements[i] = & el;
        }

    }

    target.data_length = length(source);
}

//////////////////////////////////////////////////////////////////////////////

//create Space for SkiplistElement of given height

template <typename TValue, typename TSpec>
inline SkiplistElement<TValue, TSpec> &
_skiplistAllocateElement(Map<TValue, Skiplist<TSpec> > & me,
                         unsigned char height)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef SkiplistElement<TValue, TSpec> TElement;
    typedef SkiplistNext<TValue, TSpec> TNext;

    TElement * ret;

    if (me.data_recycle[height])
    {//use recycled
        ret = me.data_recycle[height];
        me.data_recycle[height] = * reinterpret_cast<TElement **>(me.data_recycle[height]);
    }
    else
    {
        int const element_base_size = sizeof(TElement) - TSkiplist::MAX_HEIGHT * sizeof(TNext);
        int need_size = element_base_size + (height+1) * sizeof(TNext);
        int buf_size = me.data_mem_end - me.data_mem_begin;
        if (buf_size < need_size)
        {//need new memory
            if (buf_size >= (int) (element_base_size + sizeof(TNext)))
            {//link rest memory in recycle
                int rest_height = (buf_size - element_base_size) / sizeof(TNext) - 1; //must be < height, because buf_size < need_size
                * reinterpret_cast<TElement **>(me.data_mem_begin) = me.data_recycle[rest_height];
                me.data_recycle[rest_height] = reinterpret_cast<TElement *>(me.data_mem_begin);
            }
            //else: a small part of memory is wasted
            allocate(value(me.data_allocator), me.data_mem_begin, (size_t) TSkiplist::BLOCK_SIZE, TagAllocateStorage());
            me.data_mem_end = me.data_mem_begin + TSkiplist::BLOCK_SIZE;
        }
        ret = reinterpret_cast<TElement *>(me.data_mem_begin);
        me.data_mem_begin += need_size;
    }

    return *ret;
}

//////////////////////////////////////////////////////////////////////////////

//Creates New SkiplistElement
template <typename TValue, typename TSpec, typename TValue2>
inline SkiplistElement<TValue, TSpec> &
_skiplistConstructElement(Map<TValue, Skiplist<TSpec> > & me,
                          unsigned char height,
                          TValue2 const & _value)
{
    typedef SkiplistElement<TValue, TSpec> TElement;
    TElement & el = _skiplistAllocateElement(me, height);
    valueConstruct(& (el.data_value), _value);
    //no need to construct the next array

    return el;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
_skiplistDeallocateElement(Map<TValue, Skiplist<TSpec> > & me,
                           SkiplistElement<TValue, TSpec> & el,
                           unsigned char height)
{
    typedef SkiplistElement<TValue, TSpec> TElement;
    * reinterpret_cast<TElement **>(& el) = me.data_recycle[height];
    me.data_recycle[height] = reinterpret_cast<TElement *>(&el);
    //the real deallocation is done by the allocator on destruction
}

//////////////////////////////////////////////////////////////////////////////

//Destroys New SkiplistElement
template <typename TValue, typename TSpec>
inline void
_skiplistDestructElement(Map<TValue, Skiplist<TSpec> > & me,
                         SkiplistElement<TValue, TSpec> & el,
                         unsigned char height)
{
    valueDestruct(& (el.value) );
    //no need to construct the next array
    _skiplistDeallocateElement(me, el, height);
}

//////////////////////////////////////////////////////////////////////////////

//creates height for new SkiplistElement.
//increases the Skiplist height if necessary
template <typename TValue, typename TSpec>
inline unsigned char
_skiplistCreateHeight(Map<TValue, Skiplist<TSpec> > & me)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;

    unsigned char height = pickRandomNumber(me.rng, Pdf<GeometricFairCoin>());
    if (height >= TSkiplist::MAX_HEIGHT) height = TSkiplist::MAX_HEIGHT-1;

    if (height > me.data_height) me.data_height = height;

    return height;
}

template <typename TValue, typename TSpec>
inline unsigned char
_skiplistCreateHeight(Map<TValue, Skiplist<TSpec> > & me,
                      SkiplistPath<TValue, TSpec> & path) //extend path if height is increased
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;

    unsigned char height = pickRandomNumber(me.rng, Pdf<GeometricFairCoin>());
    if (height >= TSkiplist::MAX_HEIGHT) height = TSkiplist::MAX_HEIGHT-1;

    if (height > me.data_height)
    {
        for (unsigned char i = me.data_height + 1; i <= height; ++i)
        {
            path.data_elements[i] = & me.data_border;
        }
        me.data_height = height;
    }

    return height;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline unsigned char
_skiplistGetHeight(Map<TValue, Skiplist<TSpec> > & me,
                   SkiplistElement<TValue, TSpec> & el,
                   SkiplistPath<TValue, TSpec> & path)
{
    int height = me.data_height;
    for (; height > 0 ; --height)
    {
        if (path.elements[height]->data_next[height] == el) break;
    }
    return height;
}

template <typename TValue, typename TSpec>
inline unsigned char
_skiplistGetHeight(Map<TValue, Skiplist<TSpec> > & me,
                   SkiplistElement<TValue, TSpec> & el)
{
    typedef SkiplistPath<TValue, TSpec> TPath;

    TPath path;
    _skiplistFind(me, el, path);

    return _skiplistGetHeight(el, path);
}

//////////////////////////////////////////////////////////////////////////////

//Note: always store paths to the element LEFT of the actually wanted element.
//this element will be the predecessor during insertion.

//Given a key: find path to the elment that is left to:
// - the first element that has the given key, if there is any, or
// - the first element that has the smallest key larger than the given key,
//or a path to the last element, if all keys are smaller than the given key

//Given an element: find path to the element that is left to:
// - the given element, or
// - the first element that has the smallest key larger than the key of the given element,
//or a path to the last element, if all keys are smaller than the key of the given element

template <typename TValue, typename TSpec, typename TKey>
inline bool
_skiplistFindGoNext(SkiplistNext<TValue, TSpec> & next,
                    unsigned char,
                    TKey const & _key)
{
    return (key(next.data_element->data_value) < _key);
}

template <typename TValue, typename TSpec>
inline bool
_skiplistFindGoNext(SkiplistNext<TValue, TSpec> & next,
                    unsigned char,
                    SkiplistElement<TValue, TSpec> const & el)
{
    return (key(next.data_element->data_value) <= key(el.data_value)) && (next.data_element != & el);
}

template <typename TValue, typename TSpec>
inline bool
_skiplistFindGoNext(SkiplistNext<TValue, TSpec> & next,
                    unsigned char /*height*/,
                    GoEnd)
{
    return next.data_element;
//    return next.data_element->data_next[height];
}


template <typename TValue, typename TSpec, typename TFind>
inline void
_skiplistFind(Map<TValue, Skiplist<TSpec> > & me,
              TFind const & find, //can be a key or a SkiplistElement or GoEnd
              /*OUT*/ SkiplistPath<TValue, TSpec> & path)
{
    typedef SkiplistElement<TValue, TSpec> TElement;
    typedef SkiplistNext<TValue, TSpec> TNext;

    TElement * here = & me.data_border;

    for (int i = me.data_height; i >= 0; --i)
    {
        while (true)
        {
            TNext & next = here->data_next[i];
            if (!next.data_element || !_skiplistFindGoNext(next, i, find)) break;
            here = next.data_element;
        }
        path.data_elements[i] = here;
    }
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Map#find
 * @headerfile <seqan/map.h>
 * @brief Find a value in a map.
 *
 * @signature TIterator find(map, key);
 *
 * @param[in] map A map. Types: Map
 * @param[in] key A key.
 *
 * @return TIterator An iterator to the first value in <tt>map</tt> of the given key, if there is any.  An iterator
 *                   to the fist value in <tt>map</tt> with key &gt; <tt>key</tt>, otherwise.
 *
 * @see Map#value
 * @see Map#cargo
 */

template <typename TValue, typename TSpec, typename TFind>
inline typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type
find(Map<TValue, Skiplist<TSpec> > & me,
     TFind const & _find, //can be a TKey or a SkiplistElement or GoEnd
     SkiplistPath<TValue, TSpec> & path)
{
    typedef typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type TIterator;

    _skiplistFind(me, _find, path);
    return TIterator(path.data_elements[0]->data_next[0].data_element);
}
template <typename TValue, typename TSpec, typename TFind>
inline typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type
find(Map<TValue, Skiplist<TSpec> > & me,
     TFind const & _find) //can be a TKey or a SkiplistElement or GoEnd
{
    typedef SkiplistPath<TValue, TSpec> TPath;
    TPath path;
    return find(me, _find, path);
}

//////////////////////////////////////////////////////////////////////////////

//insert elements after at the position path points to
//Requirements:
// - height <= height of the skiplist
// - next must be filled at least up to height
// - el.data_next must have space for at least height many entries
template <typename TValue, typename TSpec>
inline void
_skiplistInsertElement(Map<TValue, Skiplist<TSpec> > & me,
                       SkiplistElement<TValue, TSpec> & el,
                       SkiplistPath<TValue, TSpec> & path,
                       unsigned char height)
{
    for (int i = height; i >= 0; --i)
    {
        el.data_next[i].data_element = path.data_elements[i]->data_next[i].data_element;
        path.data_elements[i]->data_next[i].data_element = & el;
    }

    ++me.data_length;
}

template <typename TValue, typename TSpec>
inline void
_skiplistInsertElement(Map<TValue, Skiplist<TSpec> > & me,
                       SkiplistElement<TValue, TSpec> & el,
                       unsigned char height)
{
    typedef SkiplistPath<TValue, TSpec> TPath;

    TPath path;
    _skiplistFind(me, key(el.data_value), path);

    _skiplistInsertElement(me, el, path, height);
}

//////////////////////////////////////////////////////////////////////////////
//creates entry if necessary

/*!
 * @fn Map#value
 * @brief Returns a value given a key.
 *
 * @signature TReference value(map, key);
 *
 * @param[in] map A map.
 * @param[in] key A key.
 *
 * @return TReference The first value in <tt>map</tt> of the given key, if there is any.  Otherwise, a new value
 *                    that is inserted to <tt>map</tt>.
 *
 * @note Do not change the key of a value in the map.
 */

template <typename TValue, typename TSpec, typename TKey2>
inline typename Value< Map<TValue, Skiplist<TSpec> > >::Type &
value(Map<TValue, Skiplist<TSpec> > & me,
      TKey2 const & _key)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef SkiplistPath<TValue, TSpec> TPath;
    typedef SkiplistElement<TValue, TSpec> TElement;
    typedef typename Iterator<TSkiplist>::Type TIterator;
    typedef typename Value<TSkiplist>::Type TValue2;

    TPath path;
    TIterator it = find(me, _key, path);
    if (it && (key(it) == _key))
    {
        return value(it);
    }
    else
    {// insert new value
        unsigned char height = _skiplistCreateHeight(me, path);
        TValue2 val_temp;
        setKey(val_temp, _key);
        TElement & el = _skiplistConstructElement(me, height, val_temp);
        _skiplistInsertElement(me, el, path, height);
        return el.data_value;
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Map#cargo
 * @brief Returns a cargo given a key.
 *
 * @signature TCargo find(map, key);
 *
 * @param[in,out] map A map.
 * @param[in]     key A key.
 *
 * @return TReturn The cargo of the first value in <tt>map</tt> of the given key, if there is any.  Otherwise, the
 *                 cargo of a new value that is inserted to <tt>map</tt>.
 */

template <typename TValue, typename TSpec, typename TKey2>
inline typename Cargo< Map<TValue, Skiplist<TSpec> > >::Type &
cargo(Map<TValue, Skiplist<TSpec> > & me,
      TKey2 const & _key)
{
    return cargo(value(me, _key));
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Map#insert
 * @brief Insert new value into map.
 *
 * @signature void insert(map, value);
 * @signature void insert(map, key, cargo);
 *
 * @param[in,out] map   A map.
 * @param[in]     value A value that is added to <tt>map</tt>.
 * @param[in]     key   A key.
 * @param[in]     cargo A cargo.
 *
 * If <tt>key</tt> and <tt>cargo</tt> are specified, a new value of that key and value is added.  If there is already a
 * value of that key in <tt>map</tt>, the value of this element is changed to <tt>cargo</tt>.
 *
 * If <tt>value</tt> is specified, and there is already a value in map of that key, than the cargo of this value is
 * changed to cargo.cargo(value).
 *
 * Use Map#add instead to insert multiple values of the same key.
 */

template <typename TValue, typename TSpec, typename TValue2>
inline void
insert(Map<TValue, Skiplist<TSpec> > & me,
       TValue2 const & _value)
{
    value(me, key(_value)) = _value;
}
template <typename TValue, typename TSpec, typename TKey2, typename TCargo2>
inline void
insert(Map<TValue, Skiplist<TSpec> > & me,
       TKey2 const & _key,
       TCargo2 const & _cargo)
{
    cargo(me, _key) = _cargo;
}

//////////////////////////////////////////////////////////////////////////////
//multiple key insert

/*!
 * @fn Map#add
 * @brief Insert another value into a multi map.
 *
 * @signature void add(map, value);
 * @signature void add(map, key, cargo);
 *
 * @param[in,out] map   A map. Types: Skiplist
 * @param[in]     value A value that is added to <tt>map</tt>.
 * @param[in]     cargo A cargo.
 * @param[in]     key   A key.
 *
 * If <tt>key</tt> and <tt>cargo</tt> are specified, a new value of that key and value is added.
 */

template <typename TValue, typename TSpec, typename TValue2>
inline void
add(Map<TValue, Skiplist<TSpec> > & me,
    TValue2 const & _value)
{
    typedef SkiplistElement<TValue, TSpec> TElement;

    unsigned char height = _skiplistCreateHeight(me);
    TElement & el = _skiplistConstructElement(me, height, _value);
    _skiplistInsertElement(me, el, height);
}
template <typename TValue, typename TSpec, typename TKey2, typename TCargo2>
inline void
add(Map<TValue, Skiplist<TSpec> > & me,
    TKey2 const & _key,
    TCargo2 const & _cargo)
{
    TValue temp_val;
    setKey(temp_val, _key);
    setCargo(temp_val, _cargo);

    add(me, temp_val);
}

//////////////////////////////////////////////////////////////////////////////

//extract element from Skiplist
template <typename TValue, typename TSpec>
inline void
_skiplistUnlinkElement(Map<TValue, Skiplist<TSpec> > & me,
                       SkiplistElement<TValue, TSpec> & el)
{
    typedef SkiplistPath<TValue, TSpec> TPath;

    TPath path;
    _skiplistFind(me, el, path);

    for (int i = me.data_height; i >= 0; --i)
    {
        if (path.data_elements[i]->data_next[i].data_element == & el)
        {
            path.data_elements[i]->data_next[i].data_element = el.data_next[i].data_element;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Map#erase
 * @brief Removes a value from a map.
 *
 * @signature void erase(map, key);
 * @signature void erase(map, iterator);
 *
 * @param[in] map      A map. Types: Map
 * @param[in] key      The key of a value in <tt>map</tt>.
 * @param[in] iterator An iterator to a value in <tt>map</tt>.
 *
 * Removes the first value in <tt>map</tt> of the given key, if there is any.
 *
 * Use @link Map#eraseAll @endlink to remove all values of the given key in a multi map.
 */

template <typename TValue, typename TSpec, typename TMap2>
inline void
erase(Map<TValue, Skiplist<TSpec> > & me,
      Iter<TMap2, SkiplistIterator> const & it)
{
    _skiplistUnlinkElement(me, * it.data_pointer);
    --me.data_length;
}

template <typename TValue, typename TSpec, typename TToRemove>
inline void
erase(Map<TValue, Skiplist<TSpec> > & me,
      TToRemove const & to_remove)
{
    typedef Map<TValue, Skiplist<TSpec> > TMap;
    typedef typename Iterator<TMap>::Type TIterator;
    TIterator it = find(me, to_remove);
    if (it && (key(it) == key(to_remove)))
    {
        erase(me, it);
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Map#eraseAll
 * @brief Removes a value from a map.
 *
 * @signature void eraseAll(map, key);
 *
 * @param[in,out] map A map. Types: Skiplist
 * @param[in]     key The key of a value in <tt>map</tt>.
 *
 * Removes all values in <tt>map</tt> of the given key, if there is any.
 */

template <typename TValue, typename TSpec, typename TToRemove>
inline void
eraseAll(Map<TValue, Skiplist<TSpec> > & me,
         TToRemove const & to_remove)
{
    typedef Map<TValue, Skiplist<TSpec> > TMap;
    typedef typename Iterator<TMap>::Type TIterator;
    TIterator it = find(me, to_remove);
    while (it && (key(it) == key(to_remove)))
    {
        TIterator it_old = it;
        ++it;
        erase(me, it_old);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(Map<TValue, Skiplist<TSpec> > & me)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;

    me.data_mem_begin = me.data_mem_end = 0;
    me.data_length = 0;
    me.data_height = 0;

    for (unsigned char i = 0; i < TSkiplist::MAX_HEIGHT; ++i)
    {
        me.data_recycle[i] = 0;
        valueConstruct(me.data_border.data_next + i, NonMinimalCtor());
    }

    clear(value(me.data_allocator));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size< Map<TValue, Skiplist<TSpec> > >::Type
length(Map<TValue, Skiplist<TSpec> > const & me)
{
    return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TIteratorSpec>
inline typename Iterator< Map<TValue, Skiplist<TSpec> >, TIteratorSpec>::Type
begin(Map<TValue, Skiplist<TSpec> > & me,
      TIteratorSpec)
{
    typedef typename Iterator< Map<TValue, Skiplist<TSpec> >, TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TValue, typename TSpec>
inline typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type
begin(Map<TValue, Skiplist<TSpec> > & me)
{
    typedef typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TIteratorSpec>
inline typename Iterator< Map<TValue, Skiplist<TSpec> >, TIteratorSpec>::Type
end(Map<TValue, Skiplist<TSpec> > &,
    TIteratorSpec)
{
    typedef typename Iterator< Map<TValue, Skiplist<TSpec> >, TIteratorSpec>::Type TIterator;
    return TIterator();
}
template <typename TValue, typename TSpec>
inline typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type
end(Map<TValue, Skiplist<TSpec> > &)
{
    typedef typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type TIterator;
    return TIterator();
}



//////////////////////////////////////////////////////////////////////////////

template <typename TCargo>
struct SkipListMapValue_
{
    template <typename TValue, typename TSpec, typename TKey2>
    static inline TCargo &
    mapValue_(Map<TValue, Skiplist<TSpec> > & me,
        TKey2 const & _key)
    {
        return cargo(me, _key);
    }
};

template <>
struct SkipListMapValue_<Nothing>
{
    template <typename TValue, typename TSpec, typename TKey2>
    static inline bool
    mapValue_(Map<TValue, Skiplist<TSpec> > & me,
        TKey2 const & _key)
    {
        return hasKey(me, _key);
    }
};

template <typename TValue, typename TSpec, typename TKey2>
inline typename MapValue< Map<TValue, Skiplist<TSpec> > >::Type
mapValue(Map<TValue, Skiplist<TSpec> > & me,
         TKey2 const & _key)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef typename Cargo<TSkiplist>::Type TCargo;
    return SkipListMapValue_<TCargo>::mapValue_(me, _key);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Map#hasKey
 * @brief Determines whether a map contains a value given key.
 *
 * @signature bool hasKey(map, key);
 *
 * @param[in] map A map. Types: Map
 * @param[in] key A key.
 *
 * @return bool <tt>true</tt>, if there is a value in <tt>map</tt> of that key, <tt>false</tt> otherwise.
 */

template <typename TValue, typename TSpec, typename TKey2>
inline bool
hasKey(Map<TValue, Skiplist<TSpec> > & me,
       TKey2 const & _key)
{
    typedef Map<TValue, Skiplist<TSpec> > TSkiplist;
    typedef typename Iterator<TSkiplist>::Type TIterator;
    TIterator it = find(me, _key);

    return (it && (key(it) == _key));
}

//////////////////////////////////////////////////////////////////////////////
// SkiplistIterator: Standard Iterator for Skiplist
//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
class Iter< TSkiplist, SkiplistIterator>
{
public:
    typedef typename SkiplistElement_<TSkiplist>::Type TElement;
    typedef typename Value<TSkiplist>::Type TValue;

    TElement * data_pointer;

    Iter()
        : data_pointer(0)
    {
    }
    Iter(Iter const & other)
        : data_pointer(other.data_pointer)
    {
    }
    Iter(TElement * ptr)
        : data_pointer(ptr)
    {
    }
    Iter(TSkiplist const & sk)
        : data_pointer(sk.data_border.data_next[0].data_element)
    {
    }
    ~Iter()
    {
    }

    Iter const &
    operator = (Iter const & other)
    {
        data_pointer = other.data_pointer;
        return *this;
    }

    // TODO(holtgrew): Why conversion to boolean?
    operator bool () const
    {
        // Using != 0 instead of implicit/explicit task here to suppress
        // performance warning C4800 in Visual Studio.
        return data_pointer != 0;
    }

    TValue *
    operator -> () const
    {
        return & data_pointer->data_value;
    }

};

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline bool
atEnd(Iter<TSkiplist, SkiplistIterator> & it)
{
    return (!it.data_pointer);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline void
goNext(Iter<TSkiplist, SkiplistIterator> & it)
{
    it.data_pointer = it.data_pointer->data_next[0].data_element;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline typename Value<TSkiplist>::Type &
value(Iter<TSkiplist, SkiplistIterator> & it)
{
    return it.data_pointer->data_value;
}
template <typename TSkiplist>
inline typename Value<TSkiplist>::Type &
value(Iter<TSkiplist, SkiplistIterator> const & it)
{
    return it.data_pointer->data_value;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline typename Key<TSkiplist>::Type const &
key(Iter<TSkiplist, SkiplistIterator> & it)
{
    return key(it.data_pointer->data_value);
}
template <typename TSkiplist>
inline typename Key<TSkiplist>::Type const &
key(Iter<TSkiplist, SkiplistIterator> const & it)
{
    return key(it.data_pointer->data_value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline typename Cargo<TSkiplist>::Type &
cargo(Iter<TSkiplist, SkiplistIterator> & it)
{
    return cargo(it.data_pointer->data_value);
}
template <typename TSkiplist>
inline typename Cargo<TSkiplist>::Type &
cargo(Iter<TSkiplist, SkiplistIterator> const & it)
{
    return cargo(it.data_pointer->data_value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline bool
operator == (Iter<TSkiplist, SkiplistIterator> const & left,
             Iter<TSkiplist, SkiplistIterator> const & right)
{
    return left.data_pointer == right.data_pointer;
}

template <typename TSkiplist>
inline bool
operator != (Iter<TSkiplist, SkiplistIterator> const & left,
             Iter<TSkiplist, SkiplistIterator> const & right)
{
    return left.data_pointer != right.data_pointer;
}

////////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

}

#endif

