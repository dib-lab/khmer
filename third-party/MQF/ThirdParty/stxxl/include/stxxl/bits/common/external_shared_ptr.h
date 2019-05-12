/***************************************************************************
 *  include/stxxl/bits/common/external_shared_ptr.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2011 Daniel Godas-Lopez <dgodas@gmail.com>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_EXTERNAL_SHARED_PTR_HEADER
#define STXXL_COMMON_EXTERNAL_SHARED_PTR_HEADER

#include <stxxl/bits/namespace.h>
#include <ostream>

STXXL_BEGIN_NAMESPACE

//! \addtogroup support
//! \{

/*!
 * This class takes a shared pointer, increments its reference count and wraps
 * it in a way that the resulting object can be copied, dumped to disk, and
 * destroyed without affecting the refcount. When the object is retrieved from
 * disk and recreated on internal memory, it will still hold a reference to the
 * same memory block and can be used right away by calling the "get" method or
 * unwrapped with the "unwrap" method to decrement the refcount.
 *
 * In the context of this template, a shared pointer is an object of a class P
 * that fulfills the following requirements:
 *
 *   - Can be copy-constructed
 *   - Has an assignment operator (so that the get method can be used)
 *   - Contains a pointer to a reference count stored outside the class
 *   - Increments the reference count on copy-construction
 *   - Decrements the reference count on destruction
 *
 * Both the Boost and c++0x implementations of shared_ptr fulfill these
 * requirements. At the moment of writing the author is not aware of any
 * implementations of shared pointers that can't be used with this wrapper.
 */
template <class P>
class external_shared_ptr
{
private:
    /*!
     * We store the pointer like this so that the refcount does not get
     * incremented when the wrapper is copy-constructed, or decremented when
     * the wrapper is destroyed.
     *
     * The whole external_shared_ptr object will be aligned by the compiler to
     * a multiple of its size. The size of the object is sizeof(P) as the
     * buffer is its only member. The buffer is placed in the class at offset 0
     * so the alignment of the stored P should be alright without any
     * additional hints.
     */
    char data[sizeof(P)];

public:
    /*!
     * This constructor needs to be defined so that the [] operator in maps and
     * hash tables works. If unwrap() or get() are called for an object
     * constructed this way the behavior is undefined.
     */
    external_shared_ptr()
    { }

    /*!
     * Copy the pointer to internal storage and increment the refcount (the
     * destructor never gets called).
     */
    external_shared_ptr(P ptr)
    {
        new (data)P(ptr);
    }

    /*!
     * Call the destructor to decrement the refcount. If this is called more
     * than once the results are undefined.
     */
    void unwrap()
    {
        P* p = reinterpret_cast<P*>((void*)data);
        p->~P();
    }

    /*!
     * If this is called after unwrap() the behaviour is undefined.
     */
    P get() const
    {
        P* p = reinterpret_cast<P*>((void*)data);
        return *p;
    }

    bool operator == (const external_shared_ptr& x) const
    {
        P* p1 = reinterpret_cast<P*>((void*)data);
        P* p2 = reinterpret_cast<P*>((void*)x.data);

        return *p1 == *p2;
    }

    //! Output contained data items
    friend std::ostream&
    operator << (std::ostream& os, const external_shared_ptr& p)
    {
        return os << p.get();
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_EXTERNAL_SHARED_PTR_HEADER
