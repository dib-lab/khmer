/***************************************************************************
 *  include/stxxl/bits/common/counting_ptr.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010-2011 Raoul Steffen <R-Steffen@gmx.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_COUNTING_PTR_HEADER
#define STXXL_COMMON_COUNTING_PTR_HEADER

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <stxxl/types>
#include <stxxl/bits/config.h>
#include <stxxl/bits/common/mutex.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup support
//! \{

/*!
 * High-performance smart pointer used as a wrapping reference counting
 * pointer.
 *
 * This smart pointer class requires two functions in the templated type: void
 * inc_reference() and void dec_reference(). These must increment and decrement
 * a reference counter inside the templated object. When initialized, the type
 * must have reference count zero. It is _not_ immediately called with
 * add_reference(). Each new object referencing the data calls add_reference()
 * and each destroying holder calls del_reference(). When the data object
 * determines that it's internal counter is zero, then it must destroy itself.
 *
 * Accompanying the counting_ptr is a const_counting_ptr and a class
 * counted_object, from which reference counted classes must be derive
 * from. The class counted_object implement all methods required for reference
 * counting.
 *
 * The whole method is more similar to boost' instrusive_ptr, but also yields
 * something resembling shared_ptr.
 */
template <class Type>
class counting_ptr
{
public:
    //! contained type.
    typedef Type element_type;

private:
    //! the pointer to the currently referenced object.
    Type* m_ptr;

protected:
    //! increment reference counter for current object.
    void inc_reference()
    { inc_reference(m_ptr); }

    //! increment reference counter of other object.
    void inc_reference(Type* o)
    { if (o) o->inc_reference(); }

    //! decrement reference counter of current object and maybe delete it.
    void dec_reference()
    { if (m_ptr && m_ptr->dec_reference()) delete m_ptr; }

public:
    //! default constructor: contains a NULL pointer.
    counting_ptr() : m_ptr(NULL)
    { }

    //! constructor with pointer: initializes new reference to ptr.
    counting_ptr(Type* ptr) : m_ptr(ptr)
    { inc_reference(); }

    //! copy-constructor: also initializes new reference to ptr.
    counting_ptr(const counting_ptr& other_ptr) : m_ptr(other_ptr)
    { inc_reference(); }

    //! assignment operator: dereference current object and acquire reference on new one.
    counting_ptr& operator = (const counting_ptr& other_ptr)
    { return operator = (other_ptr.m_ptr); }

    //! assignment to pointer: dereference current and acquire reference to new ptr.
    counting_ptr& operator = (Type* ptr)
    {
        inc_reference(ptr);
        dec_reference();
        m_ptr = ptr;
        return *this;
    }

    //! destructor: decrements reference counter in ptr.
    ~counting_ptr()
    { dec_reference(); }

    //! return the enclosed object as reference.
    Type& operator * () const
    {
        assert(m_ptr);
        return *m_ptr;
    }

    //! return the enclosed pointer.
    Type* operator -> () const
    {
        assert(m_ptr);
        return m_ptr;
    }

    //! implicit cast to the enclosed pointer.
    operator Type* () const
    { return m_ptr; }

    //! return the enclosed pointer.
    Type * get() const
    { return m_ptr; }

    //! test equality of only the pointer values.
    bool operator == (const counting_ptr& other_ptr) const
    { return m_ptr == other_ptr.m_ptr; }

    //! test inequality of only the pointer values.
    bool operator != (const counting_ptr& other_ptr) const
    { return m_ptr != other_ptr.m_ptr; }

    //! cast to bool check for a NULL pointer
    operator bool () const
    { return valid(); }

    //! test for a non-NULL pointer
    bool valid() const
    { return (m_ptr != NULL); }

    //! test for a NULL pointer
    bool empty() const
    { return (m_ptr == NULL); }

    //! if the object is referred by this counting_ptr only
    bool unique() const
    { return m_ptr && m_ptr->unique(); }

    //! make and refer a copy if the original object was shared.
    void unify()
    {
        if (m_ptr && ! m_ptr->unique())
            operator = (new Type(*m_ptr));
    }

    //! swap enclosed object with another counting pointer (no reference counts need change)
    void swap(counting_ptr& b)
    {
        std::swap(m_ptr, b.m_ptr);
    }
};

//! swap enclosed object with another counting pointer (no reference counts need change)
template <class A>
void swap(counting_ptr<A>& a1, counting_ptr<A>& a2)
{
    a1.swap(a2);
}

/*!
 * High-performance smart pointer used as a wrapping reference counting
 * pointer.
 *
 * This smart pointer class requires two functions in the templated type: void
 * inc_reference() and void dec_reference(). These must increment and decrement
 * a reference counter inside the templated object. When initialized, the type
 * must have reference count zero. It is _not_ immediately called with
 * add_reference(). Each new object referencing the data calls add_reference()
 * and each destroying holder calls del_reference(). When the data object
 * determines that it's internal counter is zero, then it must destroy itself.
 *
 * Accompanying the counting_ptr is a const_counting_ptr and a class
 * counted_object, from which reference counted classes must be derive
 * from. The class counted_object implement all methods required for reference
 * counting.
 *
 * The whole method is more similar to boost' instrusive_ptr, but also yields
 * something resembling shared_ptr.
 */
template <class Type>
class const_counting_ptr
{
public:
    //! contained type.
    typedef Type element_type;

private:
    //! the pointer to the currently referenced object.
    const Type* m_ptr;

protected:
    //! increment reference counter for current object.
    void inc_reference()
    { inc_reference(m_ptr); }

    //! increment reference counter of other object.
    void inc_reference(const Type* o)
    { if (o) o->inc_reference(); }

    //! decrement reference counter of current object and maybe delete it.
    void dec_reference()
    { if (m_ptr && m_ptr->dec_reference()) delete m_ptr; }

public:
    //! default constructor: contains a NULL pointer.
    const_counting_ptr() : m_ptr(NULL)
    { }

    //! constructor with pointer: initializes new reference to ptr.
    const_counting_ptr(const Type* ptr) : m_ptr(ptr)
    { inc_reference(); }

    //! copy-constructor: also initializes new reference to ptr.
    const_counting_ptr(const const_counting_ptr& other_ptr) : m_ptr(other_ptr)
    { inc_reference(); }

    //! constructor from non-const: also initializes new reference to ptr.
    const_counting_ptr(const counting_ptr<Type>& other_ptr) : m_ptr(other_ptr.get())
    { inc_reference(); }

    //! assignment operator: dereference current object and acquire reference on new one.
    const_counting_ptr& operator = (const const_counting_ptr& other_ptr)
    { return operator = (other_ptr.m_ptr); }

    //! assignment operator: dereference current object and acquire reference on new one.
    const_counting_ptr& operator = (const counting_ptr<Type>& other_ptr)
    { return operator = (other_ptr.get()); }

    //! assignment to pointer: dereference current and acquire reference to new ptr.
    const_counting_ptr& operator = (const Type* ptr)
    {
        inc_reference(ptr);
        dec_reference();
        m_ptr = ptr;
        return *this;
    }

    //! destructor: decrements reference counter in ptr.
    ~const_counting_ptr()
    { dec_reference(); }

    //! return the enclosed object as reference.
    const Type& operator * () const
    {
        assert(m_ptr);
        return *m_ptr;
    }

    //! return the enclosed pointer.
    const Type* operator -> () const
    {
        assert(m_ptr);
        return m_ptr;
    }

    //! implicit cast to the enclosed pointer.
    operator const Type* () const
    { return m_ptr; }

    //! return the enclosed pointer.
    const Type * get() const
    { return m_ptr; }

    //! test equality of only the pointer values.
    bool operator == (const const_counting_ptr& other_ptr) const
    { return m_ptr == other_ptr.m_ptr; }

    //! test inequality of only the pointer values.
    bool operator != (const const_counting_ptr& other_ptr) const
    { return m_ptr != other_ptr.m_ptr; }

    //! test equality of only the pointer values.
    bool operator == (const counting_ptr<Type>& other_ptr) const
    { return m_ptr == other_ptr.get(); }

    //! test inequality of only the pointer values.
    bool operator != (const counting_ptr<Type>& other_ptr) const
    { return m_ptr != other_ptr.get(); }

    //! cast to bool check for a NULL pointer
    operator bool () const
    { return m_ptr; }

    //! test for a non-NULL pointer
    bool valid() const
    { return m_ptr; }

    //! test for a NULL pointer
    bool empty() const
    { return !m_ptr; }

    //! if the object is referred by this const_counting_ptr only
    bool unique() const
    { return m_ptr && m_ptr->unique(); }

    //! swap enclosed object with another const_counting pointer (no reference counts need change)
    void swap(const_counting_ptr& b)
    {
        std::swap(m_ptr, b.m_ptr);
    }
};

//! swap enclosed object with another const_counting pointer (no reference counts need change)
template <class A>
void swap(const_counting_ptr<A>& a1, const_counting_ptr<A>& a2)
{
    a1.swap(a2);
}

/*!
 * Provides reference counting abilities for use with counting_ptr.
 *
 * Use as superclass of the actual object, this adds a reference_count
 * value. Then either use counting_ptr as pointer to manage references and
 * deletion, or just do normal new and delete.
 *
 * For thread-safe functions, use atomic_counted_object instead of this class!
 */
class counted_object
{
private:
    //! the reference count is kept mutable to all const_counting_ptr() to
    //! change the reference count.
    mutable unsigned_type m_reference_count;

public:
    //! new objects have zero reference count
    counted_object()
        : m_reference_count(0) { }

    //! coping still creates a new object with zero reference count
    counted_object(const counted_object&)
        : m_reference_count(0) { }

    //! assignment operator, leaves pointers unchanged
    counted_object& operator = (const counted_object&)
    { return *this; } // changing the contents leaves pointers unchanged

    ~counted_object()
    { assert(m_reference_count == 0); }

public:
    //! Call whenever setting a pointer to the object
    void inc_reference() const
    { ++m_reference_count; }

    //! Call whenever resetting (i.e. overwriting) a pointer to the object.
    //! IMPORTANT: In case of self-assignment, call AFTER inc_reference().
    //! \return if the object has to be deleted (i.e. if it's reference count dropped to zero)
    bool dec_reference() const
    { return (! --m_reference_count); }

    //! Test if the counted_object is referenced by only one counting_ptr.
    bool unique() const
    { return (m_reference_count == 1); }

    //! Return the number of references to this object (for debugging)
    unsigned_type get_reference_count() const
    { return m_reference_count; }
};

#if STXXL_HAVE_SYNC_ADD_AND_FETCH || STXXL_MSVC

/*!
 * Provides reference counting abilities for use with counting_ptr with atomics
 * operations.
 *
 * Use as superclass of the actual object, this adds a reference_count
 * value. Then either use counting_ptr as pointer to manage references and
 * deletion, or just do normal new and delete.
 *
 * This class does thread-safe increment and decrement using atomic operations
 * on an integral type.
 */
class atomic_counted_object
{
private:
    //! the reference count is kept mutable to all const_counting_ptr() to
    //! change the reference count.
#if STXXL_MSVC
    mutable long m_reference_count;
#else
    mutable unsigned_type m_reference_count;
#endif

public:
    //! new objects have zero reference count
    atomic_counted_object()
        : m_reference_count(0) { }

    //! coping still creates a new object with zero reference count
    atomic_counted_object(const atomic_counted_object&)
        : m_reference_count(0) { }

    //! assignment operator, leaves pointers unchanged
    atomic_counted_object& operator = (const atomic_counted_object&)
    { return *this; } // changing the contents leaves pointers unchanged

    ~atomic_counted_object()
    { assert(m_reference_count == 0); }

public:
    //! Call whenever setting a pointer to the object
    void inc_reference() const
    {
#if STXXL_MSVC
        _InterlockedIncrement(&m_reference_count);
#else
        __sync_add_and_fetch(&m_reference_count, +1);
#endif
    }

    //! Call whenever resetting (i.e. overwriting) a pointer to the object.
    //! IMPORTANT: In case of self-assignment, call AFTER inc_reference().
    //! \return if the object has to be deleted (i.e. if it's reference count dropped to zero)
    bool dec_reference() const
    {
#if STXXL_MSVC
        return (_InterlockedDecrement(&m_reference_count) == 0);
#else
        return (__sync_add_and_fetch(&m_reference_count, -1) == 0);
#endif
    }

    //! Test if the counted_object is referenced by only one counting_ptr.
    bool unique() const
    {
        return (m_reference_count == 1);
    }

    //! Return the number of references to this object (for debugging)
    unsigned_type get_reference_count() const
    {
        return m_reference_count;
    }
};

#else // no atomic intrinsics found, use mutexes (slow)

/*!
 * Provides reference counting abilities for use with counting_ptr with mutex
 * locking.
 *
 * Use as superclass of the actual object, this adds a reference_count
 * value. Then either use counting_ptr as pointer to manage references and
 * deletion, or just do normal new and delete.
 *
 * This class does thread-safe increment and decrement using scoped locks. A
 * faster version of this class is available using atomic operations.
 */
class atomic_counted_object
{
private:
    //! the reference count is kept mutable to all const_counting_ptr() to
    //! change the reference count.
    mutable unsigned_type m_reference_count;

    //! the mutex used to synchronize access to the reference counter.
    mutable mutex m_reference_count_mutex;

public:
    //! new objects have zero reference count
    atomic_counted_object()
        : m_reference_count(0) { }

    //! coping still creates a new object with zero reference count
    atomic_counted_object(const atomic_counted_object&)
        : m_reference_count(0) { }

    //! assignment operator, leaves pointers unchanged
    atomic_counted_object& operator = (const atomic_counted_object&)
    { return *this; } // changing the contents leaves pointers unchanged

    ~atomic_counted_object()
    { assert(m_reference_count == 0); }

public:
    //! Call whenever setting a pointer to the object
    void inc_reference() const
    {
        scoped_mutex_lock lock(m_reference_count_mutex);
        ++m_reference_count;
    }

    //! Call whenever resetting (i.e. overwriting) a pointer to the object.
    //! IMPORTANT: In case of self-assignment, call AFTER inc_reference().
    //! \return if the object has to be deleted (i.e. if it's reference count dropped to zero)
    bool dec_reference() const
    {
        scoped_mutex_lock lock(m_reference_count_mutex);
        return (--m_reference_count == 0);
    }

    //! Test if the counted_object is referenced by only one counting_ptr.
    bool unique() const
    {
        scoped_mutex_lock lock(m_reference_count_mutex);
        return (m_reference_count == 1);
    }

    //! Return the number of references to this object (for debugging)
    unsigned_type get_reference_count() const
    {
        scoped_mutex_lock lock(m_reference_count_mutex);
        return m_reference_count;
    }
};

#endif

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_COUNTING_PTR_HEADER
