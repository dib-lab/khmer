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

#ifndef SEQAN_HEADER_SUMLIST_MINI_H
#define SEQAN_HEADER_SUMLIST_MINI_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// MiniSumList
//////////////////////////////////////////////////////////////////////////////


template <typename TValue>
struct MiniListEntry_
{
    enum
    {
        LIMIT_0 = 1 << 6,
        LIMIT_1 = 1 << 14,
        LIMIT_2 = 1 << 30
    };

    static const unsigned char SIZES [4];

    union _union
    {
        struct struct_0
        {
            unsigned char select: 2;
            unsigned char value_0: 6;
        } _0;
        struct struct_1
        {
            unsigned short select_1: 2;
            unsigned short value_1: 14;
        } _1;
        struct struct_2
        {
            unsigned int select_2: 2;
            unsigned int value_2: 30;
        } _2;
        unsigned char value_3[sizeof(TValue)+1];
    } data;

    MiniListEntry_() {}
    MiniListEntry_(TValue val) { this->assignValue(val); }
    ~MiniListEntry_() {}

    static inline unsigned int
    sizeNeeded(TValue val)
    {
        if (val < LIMIT_1)
        {
            if (val < LIMIT_0)
            {
                return SIZES[0];
            }
            else
            {
                return SIZES[1];
            }
        }
        else
        {
            if (val < LIMIT_2)
            {
                return SIZES[2];
            }
            else
            {
                return SIZES[3];
            }
        }
    }

    inline TValue getValue() const
    {
        switch(data._0.select)
        {
        case 0: return data._0.value_0;
        case 1: return data._1.value_1;
        case 2: return data._2.value_2;
        default: return * reinterpret_cast<TValue const *>(data.value_3 + 1);
        }
    }

    inline int size() const
    {
        return SIZES[data._0.select];
    }

    inline void assignValue(TValue val)
    {
        if (val < LIMIT_1)
        {
            if (val < LIMIT_0)
            {
                data._0.select = 0;
                data._0.value_0 = val;
            }
            else
            {
                data._0.select = 1;
                data._1.value_1 = val;
            }
        }
        else
        {
            if (val < LIMIT_2)
            {
                data._0.select = 2;
                data._2.value_2 = val;
            }
            else
            {
                data._0.select = 3;
                *reinterpret_cast<TValue *>(data.value_3 + 1) = val;
            }
        }
    }
};

//____________________________________________________________________________

template <typename TValue>
const unsigned char MiniListEntry_<TValue>::SIZES [4] = {1, 2, 4, 1 + sizeof(TValue)};


//////////////////////////////////////////////////////////////////////////////

template <unsigned short SIZE = 0x0020, typename TSpec = Default>
struct MiniSumList;

struct MiniSumListValueIterator_;

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
class SumList<DIM, TValue, MiniSumList<SIZE, TSpec> >
{
public:
    typedef typename Size<SumList>::Type TSize;
    typedef SumListValues<DIM, TValue> TValues;

    unsigned char data_ [SIZE];
    TSize data_length;  //number of elements in list
    TSize data_size;  //number of bytes used in data_
    TValues data_sum; //sum of all list numbers

    SumList()
        : data_length(0)
        , data_size(0)
    {
    }

    SumList(SumList const & other)
        : data_length(other.data_length)
        , data_size(other.data_size)
        , data_sum(other.data_sum)
    {
        arrayCopyForward(other.data_, other.data_ + SIZE, data_);
    }

    ~SumList()
    {}

    SumList const &
    operator = (SumList const & other)
    {
        assign(*this, other);
        return *this;
    }


//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
struct Size< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > >
{
    typedef unsigned short Type;
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
struct Position< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > >
{
    typedef unsigned short Type;
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline void
assign(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & target,
       SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const & source)
{
    arrayCopyForward(source.data_, source.data_ + source.data_size, target.data_);
    target.data_length = source.data_length;
    target.data_size = source.data_size;
    target.data_sum = source.data_sum;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline typename Size< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > >::Type
length(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me)
{
    return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline typename Value< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > >::Type &
getSum(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me,
       unsigned int dim)
{
    return me.data_sum[dim];
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline void
clear(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me)
{
    me.data_length = 0;
    me.data_size = 0;
    clear(me.data_sum);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > >::Type
begin(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me);
}
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const>::Type
begin(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const & me)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me);
}
//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > >::Type
end(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me, GoEnd());
}
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const>::Type
end(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const& me)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me, GoEnd());
}

//////////////////////////////////////////////////////////////////////////////

//changes the value at ptr (which is at the dim-th tuple position) to new_value
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TValue2>
inline bool
_miniSumListAssignValue(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me,
                         unsigned char * ptr,
                         int dim,
                         TValue2 new_value)
{
    typedef MiniListEntry_<TValue> TEntry;

    TEntry new_entr(new_value);
    TEntry & old_entr = * reinterpret_cast<TEntry *>(ptr);

    if (new_entr.size() != old_entr.size())
    {
        int new_size = me.data_size + new_entr.size() - old_entr.size();
        if (new_size > SIZE) return false; //not enough space
        arrayCopy(ptr + old_entr.size(), me.data_ + me.data_size, ptr + new_entr.size());
    }
    me.data_sum[dim] += (new_value - old_entr.getValue());
    old_entr.assignValue(new_value);
    return true;
}


//////////////////////////////////////////////////////////////////////////////

//computes the size needed to store the values in the DIM-tupel vals
template <unsigned int DIM, typename TValue>
inline TValue
_miniSumListSizeOfValues(SumListValues<DIM, TValue> const & vals)
{
    unsigned int sum = 0;
    for (unsigned int i = 0; i < DIM; ++i)
    {
        sum += MiniListEntry_<TValue>::sizeNeeded(vals[i]);
    }
    return sum;
}


//////////////////////////////////////////////////////////////////////////////
//inserts DIM-tupel at byte_position
//returns false on capacity owerflow, true otherwise
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TPosition, typename TValue2>
inline bool
_miniSumListInsertValues(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me,
                          TPosition byte_pos,
                          TValue2 const & new_values,
                          /*OUT*/ unsigned int & new_values_size)
{
    typedef MiniListEntry_<TValue> TEntry;

    new_values_size = _miniSumListSizeOfValues<DIM>(new_values);
    unsigned int new_size = me.data_size + new_values_size;
    if (new_size > SIZE) return false; //not enough space

    if (byte_pos < me.data_size)
    { //make room
        arrayCopyBackward(me.data_ + byte_pos, me.data_ + me.data_size, me.data_ + byte_pos + new_values_size);
    }

    unsigned char * ptr = me.data_ + byte_pos;
    for (unsigned int i = 0; i < DIM; ++i)
    {
        TEntry & entr = * reinterpret_cast<TEntry *>(ptr);
        entr.assignValue(new_values[i]);
        ptr += entr.size();
    }

    me.data_size = new_size;
    me.data_sum += new_values;
    ++me.data_length;
    return true;
}


template <typename TSumList, typename TValues>
inline bool
_insertValues(TSumList & me,
              Iter<TSumList, MiniSumListValueIterator_> it,
              TValues const & new_values,
              /*OUT*/ unsigned int & new_values_size)
{
    return _miniSumListInsertValues(me, it.data_ptr - me.data_, new_values, new_values_size);
}
template <typename TSumList, typename TValue>
inline bool
_insertValues(TSumList & me,
              Iter<TSumList, MiniSumListValueIterator_> it,
              TValue const * new_values,
              /*OUT*/ unsigned int & new_values_size)
{
    SumListValues<DIMENSION<TSumList>::VALUE, typename Value<TSumList>::Type > vals(new_values);
    return _miniSumListInsertValues(me, it.data_ptr - me.data_, vals, new_values_size);
}

//////////////////////////////////////////////////////////////////////////////
//appends DIM-tupel
//returns false on capacity owerflow, true otherwise
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TValues>
inline bool
appendValues(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me,
             TValues const & new_values)
{
    unsigned int dummy;
    return _miniSumListInsertValues(me, me.data_size, new_values, dummy);
}
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TValue2>
inline bool
appendValues(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me,
             TValue2 const * new_values)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TSumList;
    SumListValues<DIMENSION<TSumList>::VALUE, typename Value<TSumList>::Type > vals(new_values);
    unsigned int dummy;
    return _miniSumListInsertValues(me, me.data_size, vals, dummy);
}

//////////////////////////////////////////////////////////////////////////////

//moves the right half of $me$ into $right$
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
inline void
splitSumList(SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & me,
             SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & right)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TSumList;
    typedef Iter<TSumList, MiniSumListValueIterator_> TIterator;
    typedef typename Size<TSumList>::Type TSize;
    typedef SumListValues<DIM, TValue> TValues;

    clear(right);

    if (me.data_size <= SIZE/2) return; //nothing to split

    TIterator it(me);
    if (atEnd(it, me)) return; //nothing to split

    TValues sum;

    TSize me_len = 0;
    while (true)
    {
        if ((it.data_ptr - me.data_) >= SIZE/2) break;

        TValues vals;
        scanValues(it, vals);
        sum += vals;
        ++me_len;
    }

    if (atEnd(it, me)) return; //nothing to split

    //split

    arrayCopyForward(it.data_ptr, me.data_ + me.data_size, right.data_);
    right.data_length = me.data_length - me_len;
    right.data_size = me.data_size - (it.data_ptr - me.data_);
    right.data_sum = me.data_sum;
    right.data_sum -= sum;

    me.data_length = me_len;
    me.data_size = it.data_ptr - me.data_;
    me.data_sum = sum;
}


//////////////////////////////////////////////////////////////////////////////
//finds tupel in which the dim-th sum would exceed val
//template <typename TSumList, typename TValue2, typename TValues>
//inline Iter<TSumList, MiniSumListValueIterator_>
//_miniSumListSearchSumList(TSumList & me,
//                           TValue2 val,
//                           int dim,
//                           /*OUT*/ TValues const & sums)
//{
//    typedef Iter<TSumList, MiniSumListValueIterator_> TIterator;
//
//    clear(sums);
//    TIterator it(me);
//
//    while (true)
//    {
//        if (atEnd(it, me)) return it;
//
//        TIterator it_next = it;
//        TValues vals;
//        scanValues(it_next, vals);
//
//        if (sums[dim] + vals[dim] >= val) return it;
//
//        it = it_next;
//    }
//}
//

//////////////////////////////////////////////////////////////////////////////
// Iterator for MiniSumList
//////////////////////////////////////////////////////////////////////////////


//struct MiniSumListIteratorRooted;
//
//template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
//class Iter<SumList<DIM, TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIteratorRooted>
//{
//public:
//    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TContainer;
//    typedef typename Size<TContainer>::Type TContainerSize;
//    typedef typename Position<TContainer>::Type TContainerPosition;
//
//    TContainer & data_container;
//    TContainerSize data_bytepos;  //position of entry in container.data_
//    TContainerPosition data_position;
//    TValue data_sum[DIM];
//
//    Iter(TContainer & cont, TContainerSize bytepos = 0, TContainerPosition pos = 0)
//        : data_container(cont)
//        , data_bytepos(bytepos)
//        , data_position(pos)
//    {
//    }
//    Iter(Iter const & other)
//        : data_container(other.data_container)
//        , data_bytepos(other.bytepos)
//        , data_position(other.data_position)
//    {
//    }
//    ~Iter()
//    {
//    }
//    inline Iter const &
//    operator = (Iter const & other)
//    {
//        data_container = other.data_container;
//        data_bytepos = other.bytepos;
//        data_position = other.data_position;
//    }
//};

//////////////////////////////////////////////////////////////////////////////

//Iterator for iterating through the single values in the minisumlist

template <typename TSumList>
class Iter<TSumList, MiniSumListValueIterator_>
{
public:
    typedef typename CopyConst_<TSumList, unsigned char>::Type TUnsignedChar;
    mutable TUnsignedChar * data_ptr;

    Iter()
        : data_ptr(0)
    {
    }
    Iter(TSumList & cont)
        : data_ptr(cont.data_)
    {
    }
    Iter(Iter const & other)
        : data_ptr(other.data_ptr)
    {
    }
    Iter(TUnsignedChar * ptr)
        : data_ptr(ptr)
    {
    }
    ~Iter()
    {
    }
    inline Iter const &
    operator = (Iter const & other)
    {
        data_ptr = other.data_ptr;
        return *this;
    }
    inline Iter const &
    operator = (TUnsignedChar * ptr)
    {
        data_ptr = ptr;
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goNext(Iter< TSumList, MiniSumListValueIterator_ > & it)
{
    typedef typename Value<TSumList>::Type TValue;
    typedef MiniListEntry_<TValue> const TEntry;
    TEntry & entr = * reinterpret_cast<TEntry *>(it.data_ptr);
    it.data_ptr += entr.size();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getValue(Iter<TSumList, MiniSumListValueIterator_ > & it)
{
    typedef typename Value<TSumList>::Type TValue;
    typedef MiniListEntry_<TValue> const TEntry;
    TEntry & entr = * reinterpret_cast<TEntry *>(it.data_ptr);
    return entr.getValue();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TValues>
inline void
scanValues(/*IN and OUT*/ Iter< TSumList, MiniSumListValueIterator_ > & it,
           /*OUT*/ TValues & values)
{
    for (unsigned int i = 0; i < DIMENSION<TSumList>::VALUE; ++i)
    {
        values[i] = getValue(it);
        goNext(it);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TSumList2>
inline bool
atEnd(Iter< TSumList, MiniSumListValueIterator_ > & it,
      TSumList2 & container)
{
    return (!it.data_ptr) || (it.data_ptr == container.data_ + container.data_size);
}

template <typename TSumList, typename TSumList2>
inline bool
atEnd(Iter< TSumList, MiniSumListValueIterator_ > & it,
      TSumList2 const & container)
{
    return (!it.data_ptr) || (it.data_ptr == container.data_ + container.data_size);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//The "real" iterator that iterates through the tuples and "knows" the sums

struct MiniSumListIterator;


template <typename TSumList>
class Iter<TSumList, MiniSumListIterator>
{
public:
    typedef Iter<TSumList, MiniSumListValueIterator_> TValueIterator;
//    typedef typename CopyConst_<TSumList, TValueIterator_>::Type TValueIterator;
    typedef typename Value<TSumList>::Type TValue;
    typedef SumListValues<DIMENSION<TSumList>::VALUE, TValue> TValues;

    TSumList * container_;
    mutable TValueIterator here_;
    mutable TValueIterator next_;
    mutable TValues values_;
    mutable TValues sums_;

    Iter()
    {
    }
    Iter(TSumList & cont)
        : container_(& cont)
        , here_(cont.data_)
        , next_(cont.data_)
    {
        if (!atEnd(next_, cont))
        {
            scanValues(next_, values_);
        }
    }
    Iter(TSumList & cont, GoEnd)
        : container_(& cont)
        , here_(cont.data_ + cont.data_size)
        , sums_(cont.data_sum)
    {
    }
    Iter(Iter const & other)
        : container_(other.container_)
        , here_(other.here_)
        , next_(other.next_)
        , values_(other.values_)
        , sums_(other.sums_)
    {
    }

    ~Iter()
    {
    }
    inline Iter const &
    operator = (Iter const & other)
    {
        container_ = other.container_;
        here_ = other.here_;
        next_ = other.next_;
        values_ = other.values_;
        sums_ = other.sums_;
        return *this;
    }
};


//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TIteratorSpec>
struct Iterator< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> >, TIteratorSpec>
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TSumList_;
    typedef Iter<TSumList_, MiniSumListIterator> Type;
};
template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TIteratorSpec>
struct Iterator< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const, TIteratorSpec>
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > const TSumList_;
    typedef Iter<TSumList_, MiniSumListIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

//template <typename TSumList>
//struct Values< Iter<TSumList, MiniSumListIterator> >
//{
//    typedef typename Value<TSumList>::Type TValue;
//    typedef SumListValues<DIMENSION<TSumList>::VALUE, TValue> Type;
//};

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
clear(Iter<TSumList, MiniSumListIterator> & it)
{
    it.here_ = it.next_ = 0;
    it.container_ = 0;
    clear(it.values_);
    clear(it.sums_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goNext(Iter<TSumList, MiniSumListIterator> & it)
{
    it.sums_ += it.values_;
    it.here_ = it.next_;
    if (!atEnd(it.next_, *it.container_))
    {
        scanValues(it.next_, it.values_);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goBegin(Iter<TSumList, MiniSumListIterator> & it)
{
    clear(it.sums_);
    it.here_ = it.next_ = it.container_->data_;
    scanValues(it.next_, it.values_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goBeforeEnd(Iter<TSumList, MiniSumListIterator> & it)
{
    goBegin(it);
    while (!atEnd(it.next_, *it.container_))
    {
        goNext(it);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline bool
atEnd(Iter<TSumList, MiniSumListIterator> & it)
{
    return atEnd(it.here_, * it.container_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getValue(Iter<TSumList, MiniSumListIterator > & it,
         int dim)
{
    return it.values_[dim];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getSum(Iter<TSumList, MiniSumListIterator > & it,
         int dim)
{
    return it.sums_[dim];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TValue>
inline void
searchSumList(Iter< TSumList, MiniSumListIterator > & it,
              TValue const & val,
              int dim)
{
    goBegin(it);
    while (!atEnd(it) && ((getSum(it, dim) + getValue(it, dim))<= val)) // SEARCH SEMANTICS
    {
        goNext(it);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec, typename TValue2>
inline bool
assignValue(Iter<SumList<DIM, TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it,
            int dim,
            TValue2 val)
{
    typedef SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > TSumList;
    typedef Iter<TSumList, MiniSumListValueIterator_> TValueIterator;
    typedef MiniListEntry_<TValue> TEntry;

    //find a pointer to the value
    TValueIterator vit = it.here_;
    for (int i = 0; i < dim; ++i)
    {
        goNext(vit);
    }

    //save old value
    TValue old_val = getValue(vit);

    //create TEntry accessors
    TEntry new_entr(val);
    TEntry & old_entr = * reinterpret_cast<TEntry *>(vit.data_ptr);

    int diff = new_entr.size() -  old_entr.size();
    if (diff)
    {//entries take different spaces:

        int new_size = it.container_->data_size + diff;
        if (new_size > SIZE) return false; //not enough space

        //make room!
        arrayCopy(vit.data_ptr + old_entr.size(), it.container_->data_ + it.container_->data_size, vit.data_ptr + new_entr.size());
        it.container_->data_size += diff;

        //adjust next_ pointer
        it.next_.data_ptr += diff;
    }
    //adjust sum
    it.container_->data_sum[dim] += (val - old_val);

    //save new value
    old_entr.assignValue(val);
    it.values_[dim] = val;

    return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TValues>
inline bool
insertValues(Iter<TSumList, MiniSumListIterator > & it,
             TValues const & vals)
{
    unsigned int new_values_size;
    bool ret = _insertValues(*it.container_, it.here_, vals, new_values_size);
    if (ret)
    {
        it.next_.data_ptr = it.here_.data_ptr + new_values_size;
        it.values_ = vals;
    }
    return ret;
}
template <typename TSumList, typename TValue>
inline bool
insertValues(Iter<TSumList, MiniSumListIterator > & it,
             TValue const * p_vals)
{
    SumListValues<DIMENSION<TSumList>::VALUE, typename Value<TSumList>::Type > vals(p_vals);
    return insertValues(it, vals);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
removeValues(Iter<TSumList, MiniSumListIterator > & it)
{
    //adjust sums
    for (int i = 0; i < DIMENSION<TSumList>::VALUE; ++i)
    {
        it.container_->data_sum[i] -= it.values_[i];
    }
    if (!atEnd(it.next_, * it.container_))
    {//move the rest
        arrayCopyForward(it.next_.data_ptr, it.container_->data_ + it.container_->data_size, it.here_.data_ptr);
    }
    //adjust size and length
    it.container_->data_size -= (it.next_.data_ptr - it.here_.data_ptr);
    --it.container_->data_length;

    //find new next pointer
    it.next_ = it.here_;
    if (!atEnd(it.next_, *it.container_))
    {
        scanValues(it.next_, it.values_);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TSumList2>
inline bool
operator == (Iter<TSumList, MiniSumListIterator> const & left,
             Iter<TSumList2, MiniSumListIterator> const & right)
{
    return left.here_.data_ptr == right.here_.data_ptr;
}

template <typename TSumList, typename TSumList2>
inline bool
operator != (Iter<TSumList, MiniSumListIterator> const & left,
             Iter<TSumList2, MiniSumListIterator> const & right)
{
    return left.here_.data_ptr != right.here_.data_ptr;
}

//////////////////////////////////////////////////////////////////////////////

//assignValues

//////////////////////////////////////////////////////////////////////////////

//template <unsigned int DIM, typename TValue, unsigned short SIZE, typename TSpec>
//inline bool
//atEnd(Iter< SumList<DIM, TValue, MiniSumList<SIZE, TSpec> >, MiniSumListValueIterator_ > & it,
//      SumList<DIM, TValue, MiniSumList<SIZE, TSpec> > & container)
//{
//    return it.data_ptr == container.data_ + container.data_size;
//}


//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
