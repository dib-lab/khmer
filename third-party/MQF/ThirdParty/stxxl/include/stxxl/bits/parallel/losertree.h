/***************************************************************************
 *  include/stxxl/bits/parallel/losertree.h
 *
 *  Many generic loser tree variants.
 *  Extracted from MCSTL - http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014-2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_LOSERTREE_HEADER
#define STXXL_PARALLEL_LOSERTREE_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/parallel/base.h>
#include <functional>

STXXL_BEGIN_NAMESPACE

namespace parallel {

/**
 * Guarded loser tree/tournament tree, either copying the whole element into
 * the tree structure, or looking up the element via the index.
 *
 * This is a base class for the LoserTreeCopy\<true> and \<false> classes.
 *
 * Guarding is done explicitly through one flag sup per element, inf is not
 * needed due to a better initialization routine.  This is a well-performing
 * variant.
 *
 * \tparam ValueType the element type
 * \tparam Comparator comparator to use for binary comparisons.
 */
template <typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreeCopyBase
{
public:
    //! size of counters and array indexes
    typedef unsigned int size_type;
    //! type of the source field
    typedef int source_type;

protected:
    //! Internal representation of a loser tree player/node
    struct Loser
    {
        //! flag, true iff is a virtual maximum sentinel
        bool sup;
        //! index of source
        source_type source;
        //! copy of key value of the element in this node
        ValueType key;
    };

    //! number of nodes
    const size_type ik;
    //! log_2(ik) next greater power of 2
    const size_type k;
    //! array containing loser tree nodes
    Loser* losers;
    //! the comparator object
    Comparator comp;
    //! still have to construct keys
    bool first_insert;

public:
    LoserTreeCopyBase(size_type _k,
                      Comparator _comp = std::less<ValueType>())
        : ik(_k),
          k(round_up_to_power_of_two(ik)),
          comp(_comp),
          first_insert(true)
    {
        // avoid default-constructing losers[].key
        losers = static_cast<Loser*>(operator new (2 * k * sizeof(Loser)));

        for (size_type i = ik - 1; i < k; ++i)
        {
            losers[i + k].sup = true;
            losers[i + k].source = (source_type)(-1);
        }
    }

    ~LoserTreeCopyBase()
    {
        for (size_type i = 0; i < (2 * k); ++i)
            losers[i].~Loser();

        delete losers;
    }

    void print(std::ostream& os)
    {
        for (size_type i = 0; i < (k * 2); i++)
            os << i << "    " << losers[i].key << " from " << losers[i].source << ",  " << losers[i].sup << "\n";
    }

    //! return the index of the player with the smallest element.
    source_type get_min_source()
    {
        return losers[0].source;
    }

    /**
     * Initializes the player source with the element key.
     *
     * \param key the element to insert
     * \param source index of the player
     * \param sup flag that determines whether the value to insert is an
     *   explicit supremum sentinel.
     */
    void insert_start(const ValueType& key, source_type source, bool sup)
    {
        size_type pos = k + source;

        losers[pos].sup = sup;
        losers[pos].source = source;

        if (UNLIKELY(first_insert))
        {
            // copy construct all keys from this first key
            for (size_type i = 0; i < (2 * k); ++i)
                new (&(losers[i].key))ValueType(key);
            first_insert = false;
        }
        else
            losers[pos].key = key;
    }

    /**
     * Computes the winner of the competition at player root.  Called
     * recursively (starting at 0) to build the initial tree.
     *
     * \param root index of the game to start.
     */
    size_type init_winner(size_type root)
    {
        if (root >= k)
        {
            return root;
        }
        else
        {
            size_type left = init_winner(2 * root);
            size_type right = init_winner(2 * root + 1);
            if (losers[right].sup ||
                (!losers[left].sup && !comp(losers[right].key, losers[left].key)))
            {                   //left one is less or equal
                losers[root] = losers[right];
                return left;
            }
            else
            {                   //right one is less
                losers[root] = losers[left];
                return right;
            }
        }
    }

    void init()
    {
        losers[0] = losers[init_winner(1)];
    }
};

/**
 * Guarded loser tree/tournament tree, either copying the whole element into
 * the tree structure, or looking up the element via the index.
 *
 * Unstable specialization of LoserTreeCopyBase.
 *
 * Guarding is done explicitly through one flag sup per element, inf is not
 * needed due to a better initialization routine.  This is a well-performing
 * variant.
 *
 * \tparam ValueType the element type
 * \tparam Comparator comparator to use for binary comparisons.
 */
template <bool Stable /* == false */,
          typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreeCopy : public LoserTreeCopyBase<ValueType, Comparator>
{
public:
    typedef LoserTreeCopyBase<ValueType, Comparator> base_type;

    typedef typename base_type::size_type size_type;
    typedef typename base_type::source_type source_type;
    using base_type::k;
    using base_type::losers;
    using base_type::comp;

    LoserTreeCopy(size_type _k, Comparator _comp = std::less<ValueType>())
        : base_type(_k, _comp)
    { }

    // do not pass const reference since key will be used as local variable
    void delete_min_insert(ValueType key, bool sup)
    {
        using std::swap;

        source_type source = losers[0].source;
        for (size_type pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            // the smaller one gets promoted
            if (sup ||
                (!losers[pos].sup && comp(losers[pos].key, key)))
            {
                // the other one is smaller
                swap(losers[pos].sup, sup);
                swap(losers[pos].source, source);
                swap(losers[pos].key, key);
            }
        }

        losers[0].sup = sup;
        losers[0].source = source;
        losers[0].key = key;
    }
};

/**
 * Guarded loser tree/tournament tree, either copying the whole element into
 * the tree structure, or looking up the element via the index.
 *
 * Stable specialization of LoserTreeCopyBase.
 *
 * Guarding is done explicitly through one flag sup per element, inf is not
 * needed due to a better initialization routine.  This is a well-performing
 * variant.
 *
 * \tparam ValueType the element type
 * \tparam Comparator comparator to use for binary comparisons.
 */
template <typename ValueType, typename Comparator>
class LoserTreeCopy</* Stable == */ true, ValueType, Comparator>
    : public LoserTreeCopyBase<ValueType, Comparator>
{
public:
    typedef LoserTreeCopyBase<ValueType, Comparator> base_type;

    typedef typename base_type::size_type size_type;
    typedef typename base_type::source_type source_type;
    using base_type::k;
    using base_type::losers;
    using base_type::comp;

    LoserTreeCopy(size_type _k, Comparator _comp = std::less<ValueType>())
        : base_type(_k, _comp)
    { }

    // do not pass const reference since key will be used as local variable
    void delete_min_insert(ValueType key, bool sup)
    {
        using std::swap;

        source_type source = losers[0].source;
        for (size_type pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            if ((sup && (!losers[pos].sup || losers[pos].source < source)) ||
                (!sup && !losers[pos].sup &&
                 ((comp(losers[pos].key, key)) ||
                  (!comp(key, losers[pos].key) && losers[pos].source < source))))
            {
                // the other one is smaller
                swap(losers[pos].sup, sup);
                swap(losers[pos].source, source);
                swap(losers[pos].key, key);
            }
        }

        losers[0].sup = sup;
        losers[0].source = source;
        losers[0].key = key;
    }
};

/** Guarded loser tree, either copying the whole element into the tree structure, or looking up the element via the index.
 *
 *  Guarding is done explicitly through one flag sup per element, inf is not needed due to a better initialization routine.
 *  This is a well-performing variant.
 */
template <typename T, typename Comparator = std::less<T> >
class LoserTreeReference
{
#undef COPY
#ifdef COPY
        #define KEY(i) losers[i].key
        #define KEY_SOURCE(i) key
#else
        #define KEY(i) keys[losers[i].source]
        #define KEY_SOURCE(i) keys[i]
#endif

private:
    struct Loser
    {
        bool sup;
        int source;
#ifdef COPY
        T key;
#endif
    };

    unsigned int ik, k;
    Loser* losers;
#ifndef COPY
    T* keys;
#endif
    Comparator comp;

public:
    LoserTreeReference(unsigned int _k, Comparator _comp = std::less<T>()) : comp(_comp)
    {
        ik = _k;
        k = round_up_to_power_of_two(ik);
        losers = new Loser[k * 2];
#ifndef COPY
        keys = new T[ik];
#endif
        for (unsigned int i = ik - 1; i < k; i++)
            losers[i + k].sup = true;
    }

    ~LoserTreeReference()
    {
        delete[] losers;
#ifndef COPY
        delete[] keys;
#endif
    }

    void print(std::ostream& os)
    {
        for (unsigned int i = 0; i < (k * 2); i++)
            os << i << "    " << KEY(i) << " from " << losers[i].source << ",  " << losers[i].sup << "\n";
    }

    int get_min_source()
    {
        return losers[0].source;
    }

    void insert_start(T key, int source, bool sup)
    {
        unsigned int pos = k + source;

        losers[pos].sup = sup;
        losers[pos].source = source;
        KEY(pos) = key;
    }

    unsigned int init_winner(unsigned int root)
    {
        if (root >= k)
        {
            return root;
        }
        else
        {
            unsigned int left = init_winner(2 * root);
            unsigned int right = init_winner(2 * root + 1);
            if (losers[right].sup ||
                (!losers[left].sup && !comp(KEY(right), KEY(left))))
            {                   //left one is less or equal
                losers[root] = losers[right];
                return left;
            }
            else
            {                   //right one is less
                losers[root] = losers[left];
                return right;
            }
        }
    }

    void init()
    {
        losers[0] = losers[init_winner(1)];
    }

    void delete_min_insert(T /* key */, bool sup)
    {
        using std::swap;

        int source = losers[0].source;
        for (unsigned int pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            //the smaller one gets promoted
            if (sup ||
                (!losers[pos].sup && comp(KEY(pos), KEY_SOURCE(source))))
            {                   //the other one is smaller
                swap(losers[pos].sup, sup);
                swap(losers[pos].source, source);
#ifdef COPY
                swap(KEY(pos), KEY_SOURCE(source));
#endif
            }
        }

        losers[0].sup = sup;
        losers[0].source = source;
#ifdef COPY
        KEY(0) = KEY_SOURCE(source);
#endif
    }

    void insert_start_stable(T key, int source, bool sup)
    {
        return insert_start(key, source, sup);
    }

    unsigned int init_winner_stable(unsigned int root)
    {
        if (root >= k)
        {
            return root;
        }
        else
        {
            unsigned int left = init_winner(2 * root);
            unsigned int right = init_winner(2 * root + 1);
            if (losers[right].sup ||
                (!losers[left].sup && !comp(KEY(right), KEY(left))))
            {                   //left one is less or equal
                losers[root] = losers[right];
                return left;
            }
            else
            {                   //right one is less
                losers[root] = losers[left];
                return right;
            }
        }
    }

    void init_stable()
    {
        losers[0] = losers[init_winner_stable(1)];
    }

    void delete_min_insert_stable(T /* key */, bool sup)
    {
        using std::swap;

        int source = losers[0].source;
        for (unsigned int pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            //the smaller one gets promoted, ties are broken by source
            if ((sup && (!losers[pos].sup || losers[pos].source < source)) ||
                (!sup && !losers[pos].sup &&
                 ((comp(KEY(pos), KEY_SOURCE(source))) ||
                  (!comp(KEY_SOURCE(source), KEY(pos)) && losers[pos].source < source))))
            {                   //the other one is smaller
                swap(losers[pos].sup, sup);
                swap(losers[pos].source, source);
#ifdef COPY
                swap(KEY(pos), KEY_SOURCE(source));
#endif
            }
        }

        losers[0].sup = sup;
        losers[0].source = source;
#ifdef COPY
        KEY(0) = KEY_SOURCE(source);
#endif
    }
};
#undef KEY
#undef KEY_SOURCE

/**
 * Guarded loser tree, using pointers to the elements instead of copying them
 * into the tree nodes.
 *
 * This is a base class for the LoserTreePointer\<true> and \<false> classes.
 *
 * Guarding is done explicitly through one flag sup per element, inf is not
 * needed due to a better initialization routine.  This is a well-performing
 * variant.
 */
template <typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreePointerBase
{
public:
    //! size of counters and array indexes
    typedef typename LoserTreeCopyBase<ValueType, Comparator>
        ::size_type size_type;
    //! type of the source field
    typedef typename LoserTreeCopyBase<ValueType, Comparator>
        ::source_type source_type;

protected:
    //! Internal representation of a loser tree player/node
    struct Loser
    {
        //! flag, true iff is a virtual maximum sentinel
        bool sup;
        //! index of source
        source_type source;
        //! pointer to key value of the element in this node
        const ValueType* keyp;
    };

    //! number of nodes
    const size_type ik;
    //! log_2(ik) next greater power of 2
    const size_type k;
    //! array containing loser tree nodes
    Loser* losers;
    //! the comparator object
    Comparator comp;

public:
    LoserTreePointerBase(size_type _k,
                         Comparator _comp = std::less<ValueType>())
        : ik(_k),
          k(round_up_to_power_of_two(ik)),
          losers(new Loser[k * 2]),
          comp(_comp)
    {
        for (size_type i = ik - 1; i < k; i++)
        {
            losers[i + k].sup = true;
            losers[i + k].source = (source_type)(-1);
        }
    }

    ~LoserTreePointerBase()
    {
        delete[] losers;
    }

    void print(std::ostream& os)
    {
        for (size_type i = 0; i < (k * 2); i++)
            os << i << "    " << losers[i].keyp << " from " << losers[i].source << ",  " << losers[i].sup << "\n";
    }

    //! return the index of the player with the smallest element.
    source_type get_min_source()
    {
        return losers[0].source;
    }

    /**
     * Initializes the player source with the element key.
     *
     * \param key the element to insert
     * \param source index of the player
     * \param sup flag that determines whether the value to insert is an
     *   explicit supremum sentinel.
     */
    void insert_start(const ValueType& key, source_type source, bool sup)
    {
        size_type pos = k + source;

        losers[pos].sup = sup;
        losers[pos].source = source;
        losers[pos].keyp = &key;
    }

    /**
     * Computes the winner of the competition at player root.  Called
     * recursively (starting at 0) to build the initial tree.
     *
     * \param root index of the game to start.
     */
    size_type init_winner(size_type root)
    {
        if (root >= k)
        {
            return root;
        }
        else
        {
            size_type left = init_winner(2 * root);
            size_type right = init_winner(2 * root + 1);
            if (losers[right].sup ||
                (!losers[left].sup && !comp(*losers[right].keyp, *losers[left].keyp)))
            {                   //left one is less or equal
                losers[root] = losers[right];
                return left;
            }
            else
            {                   //right one is less
                losers[root] = losers[left];
                return right;
            }
        }
    }

    void init()
    {
        losers[0] = losers[init_winner(1)];
    }
};

/**
 * Guarded loser tree, using pointers to the elements instead of copying them
 * into the tree nodes.
 *
 * Unstable specialization of LoserTreeCopyBase.
 *
 * Guarding is done explicitly through one flag sup per element, inf is not
 * needed due to a better initialization routine.  This is a well-performing
 * variant.
 */
template <bool Stable /* == false */,
          typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreePointer : public LoserTreePointerBase<ValueType, Comparator>
{
public:
    typedef LoserTreePointerBase<ValueType, Comparator> base_type;

    typedef typename base_type::size_type size_type;
    typedef typename base_type::source_type source_type;
    using base_type::k;
    using base_type::losers;
    using base_type::comp;

    LoserTreePointer(size_type _k, Comparator _comp = std::less<ValueType>())
        : base_type(_k, _comp)
    { }

    void delete_min_insert(const ValueType& key, bool sup)
    {
        using std::swap;

        const ValueType* keyp = &key;
        source_type source = losers[0].source;
        for (size_type pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            //the smaller one gets promoted
            if (sup ||
                (!losers[pos].sup && comp(*losers[pos].keyp, *keyp)))
            {                   //the other one is smaller
                swap(losers[pos].sup, sup);
                swap(losers[pos].source, source);
                swap(losers[pos].keyp, keyp);
            }
        }

        losers[0].sup = sup;
        losers[0].source = source;
        losers[0].keyp = keyp;
    }
};

/**
 * Guarded loser tree, using pointers to the elements instead of copying them
 * into the tree nodes.
 *
 * Unstable specialization of LoserTreeCopyBase.
 *
 * Guarding is done explicitly through one flag sup per element, inf is not
 * needed due to a better initialization routine.  This is a well-performing
 * variant.
 */
template <typename ValueType, typename Comparator>
class LoserTreePointer</* Stable == */ true, ValueType, Comparator>
    : public LoserTreePointerBase<ValueType, Comparator>
{
public:
    typedef LoserTreePointerBase<ValueType, Comparator> base_type;

    typedef typename base_type::size_type size_type;
    typedef typename base_type::source_type source_type;
    using base_type::k;
    using base_type::losers;
    using base_type::comp;

    LoserTreePointer(size_type _k, Comparator _comp = std::less<ValueType>())
        : base_type(_k, _comp)
    { }

    void delete_min_insert(const ValueType& key, bool sup)
    {
        using std::swap;

        const ValueType* keyp = &key;
        source_type source = losers[0].source;
        for (size_type pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            //the smaller one gets promoted, ties are broken by source
            if ((sup && (!losers[pos].sup || losers[pos].source < source)) ||
                (!sup && !losers[pos].sup &&
                 ((comp(*losers[pos].keyp, *keyp)) ||
                  (!comp(*keyp, *losers[pos].keyp) && losers[pos].source < source))))
            {                   //the other one is smaller
                swap(losers[pos].sup, sup);
                swap(losers[pos].source, source);
                swap(losers[pos].keyp, keyp);
            }
        }

        losers[0].sup = sup;
        losers[0].source = source;
        losers[0].keyp = keyp;
    }
};

/**
 * Unguarded loser tree, copying the whole element into the tree structure.
 *
 * This is a base class for the LoserTreeCopyUnguarded\<true> and \<false>
 * classes.
 *
 * No guarding is done, therefore not a single input sequence must run empty.
 * This is a very fast variant.
 */
template <typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreeCopyUnguardedBase : private noncopyable
{
protected:
    //! Internal representation of a loser tree player/node
    struct Loser
    {
        //! index of source
        int source;
        //! copy of key value of the element in this node
        ValueType key;
    };

    //! number of nodes
    unsigned int ik;
    //! log_2(ik) next greater power of 2
    unsigned int k;
    //! array containing loser tree nodes
    Loser* losers;
    //! the comparator object
    Comparator comp;

public:
    LoserTreeCopyUnguardedBase(unsigned int _k, const ValueType& _sentinel,
                               Comparator _comp = std::less<ValueType>())
        : ik(_k),
          k(round_up_to_power_of_two(ik)),
          losers(new Loser[k * 2]),
          comp(_comp)
    {
        for (unsigned int i = 0; i < 2 * k; i++)
        {
            losers[i].source = -1;
            losers[i].key = _sentinel;
        }
    }

    ~LoserTreeCopyUnguardedBase()
    {
        delete[] losers;
    }

    void print(std::ostream& os)
    {
        for (unsigned int i = 0; i < k + ik; i++)
            os << i << "    " << losers[i].key << " from " << losers[i].source << "\n";
    }

    //! return the index of the player with the smallest element.
    int get_min_source()
    {
        assert(losers[0].source != -1 && "Data underrun in unguarded merging.");
        return losers[0].source;
    }

    void insert_start(const ValueType& key, int source)
    {
        unsigned int pos = k + source;

        losers[pos].source = source;
        losers[pos].key = key;
    }

    unsigned int init_winner(unsigned int root)
    {
        if (root >= k)
        {
            return root;
        }
        else
        {
            unsigned int left = init_winner(2 * root);
            unsigned int right = init_winner(2 * root + 1);
            if (!comp(losers[right].key, losers[left].key))
            {                   //left one is less or equal
                losers[root] = losers[right];
                return left;
            }
            else
            {                   //right one is less
                losers[root] = losers[left];
                return right;
            }
        }
    }

    void init()
    {
        losers[0] = losers[init_winner(1)];
    }
};

template <bool Stable /* == false */,
          typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreeCopyUnguarded : public LoserTreeCopyUnguardedBase<ValueType, Comparator>
{
protected:
    typedef LoserTreeCopyUnguardedBase<ValueType, Comparator> base_type;

    using base_type::k;
    using base_type::losers;
    using base_type::comp;

public:
    LoserTreeCopyUnguarded(unsigned int _k, const ValueType& _sentinel,
                           Comparator _comp = std::less<ValueType>())
        : base_type(_k, _sentinel, _comp)
    { }

    // do not pass const reference since key will be used as local variable
    void delete_min_insert(ValueType key)
    {
        using std::swap;

        int source = losers[0].source;
        for (unsigned int pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            // the smaller one gets promoted
            if (comp(losers[pos].key, key))
            {
                // the other one is smaller
                swap(losers[pos].source, source);
                swap(losers[pos].key, key);
            }
        }

        losers[0].source = source;
        losers[0].key = key;
    }
};

template <typename ValueType, typename Comparator>
class LoserTreeCopyUnguarded</* Stable == */ true, ValueType, Comparator>
    : public LoserTreeCopyUnguardedBase<ValueType, Comparator>
{
protected:
    typedef LoserTreeCopyUnguardedBase<ValueType, Comparator> base_type;

    using base_type::k;
    using base_type::losers;
    using base_type::comp;

public:
    LoserTreeCopyUnguarded(unsigned int _k, const ValueType& _sentinel,
                           Comparator _comp = std::less<ValueType>())
        : base_type(_k, _sentinel, _comp)
    { }

    // do not pass const reference since key will be used as local variable
    void delete_min_insert(ValueType key)
    {
        using std::swap;

        int source = losers[0].source;
        for (unsigned int pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            if (!comp(key, losers[pos].key) && losers[pos].source < source)
            {
                // the other one is smaller
                swap(losers[pos].source, source);
                swap(losers[pos].key, key);
            }
        }

        losers[0].source = source;
        losers[0].key = key;
    }
};

/**
 * Unguarded loser tree, keeping only pointers to the elements in the tree
 * structure.
 *
 * This is a base class for the LoserTreePointerUnguarded\<true> and \<false>
 * classes.
 *
 * No guarding is done, therefore not a single input sequence must run empty.
 * This is a very fast variant.
 */
template <typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreePointerUnguardedBase : private noncopyable
{
protected:
    //! Internal representation of a loser tree player/node
    struct Loser
    {
        //! index of source
        int source;
        //! copy of key value of the element in this node
        const ValueType* keyp;
    };

    //! number of nodes
    unsigned int ik;
    //! log_2(ik) next greater power of 2
    unsigned int k;
    //! array containing loser tree nodes
    Loser* losers;
    //! the comparator object
    Comparator comp;

public:
    LoserTreePointerUnguardedBase(unsigned int _k, const ValueType& _sentinel,
                                  Comparator _comp = std::less<ValueType>())
        : ik(_k),
          k(round_up_to_power_of_two(ik)),
          losers(new Loser[k * 2]),
          comp(_comp)
    {
        for (unsigned int i = ik - 1; i < k; i++)
        {
            losers[i + k].source = -1;
            losers[i + k].keyp = &_sentinel;
        }
    }

    ~LoserTreePointerUnguardedBase()
    {
        delete[] losers;
    }

    void print(std::ostream& os)
    {
        for (unsigned int i = 0; i < k + ik; i++)
            os << i << "    " << *losers[i].keyp << " from " << losers[i].source << "\n";
    }

    int get_min_source()
    {
        return losers[0].source;
    }

    void insert_start(const ValueType& key, int source)
    {
        unsigned int pos = k + source;

        losers[pos].source = source;
        losers[pos].keyp = &key;
    }

    unsigned int init_winner(unsigned int root)
    {
        if (root >= k)
        {
            return root;
        }
        else
        {
            unsigned int left = init_winner(2 * root);
            unsigned int right = init_winner(2 * root + 1);
            if (!comp(*losers[right].keyp, *losers[left].keyp))
            {                   //left one is less or equal
                losers[root] = losers[right];
                return left;
            }
            else
            {                   //right one is less
                losers[root] = losers[left];
                return right;
            }
        }
    }

    void init()
    {
        losers[0] = losers[init_winner(1)];
    }
};

template <bool Stable /* == false */,
          typename ValueType, typename Comparator = std::less<ValueType> >
class LoserTreePointerUnguarded
    : public LoserTreePointerUnguardedBase<ValueType, Comparator>
{
protected:
    typedef LoserTreePointerUnguardedBase<ValueType, Comparator> base_type;

    using base_type::k;
    using base_type::losers;
    using base_type::comp;

public:
    LoserTreePointerUnguarded(unsigned int _k, const ValueType& _sentinel,
                              Comparator _comp = std::less<ValueType>())
        : base_type(_k, _sentinel, _comp)
    { }

    void delete_min_insert(const ValueType& key)
    {
        using std::swap;

        const ValueType* keyp = &key;
        int source = losers[0].source;
        for (unsigned int pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            //the smaller one gets promoted
            if (comp(*losers[pos].keyp, *keyp))
            {
                //the other one is smaller
                swap(losers[pos].source, source);
                swap(losers[pos].keyp, keyp);
            }
        }

        losers[0].source = source;
        losers[0].keyp = keyp;
    }
};

template <typename ValueType, typename Comparator>
class LoserTreePointerUnguarded</* Stable == */ true, ValueType, Comparator>
    : public LoserTreePointerUnguardedBase<ValueType, Comparator>
{
protected:
    typedef LoserTreePointerUnguardedBase<ValueType, Comparator> base_type;

    using base_type::k;
    using base_type::losers;
    using base_type::comp;

public:
    LoserTreePointerUnguarded(unsigned int _k, const ValueType& _sentinel,
                              Comparator _comp = std::less<ValueType>())
        : base_type(_k, _sentinel, _comp)
    { }

    void delete_min_insert(const ValueType& key)
    {
        using std::swap;

        const ValueType* keyp = &key;
        int source = losers[0].source;
        for (unsigned int pos = (k + source) / 2; pos > 0; pos /= 2)
        {
            //the smaller one gets promoted, ties are broken by source
            if (comp(*losers[pos].keyp, *keyp) ||
                (!comp(*keyp, *losers[pos].keyp) && losers[pos].source < source))
            {
                //the other one is smaller
                swap(losers[pos].source, source);
                swap(losers[pos].keyp, keyp);
            }
        }

        losers[0].source = source;
        losers[0].keyp = keyp;
    }
};

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_LOSERTREE_HEADER
