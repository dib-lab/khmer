/***************************************************************************
 *  include/stxxl/bits/common/new_alloc.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2006 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007, 2008, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_NEW_ALLOC_HEADER
#define STXXL_COMMON_NEW_ALLOC_HEADER

#include <memory>
#include <limits>
#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

template <class Type>
class new_alloc;

template <typename Type, typename Rebind>
struct new_alloc_rebind;

template <typename Type>
struct new_alloc_rebind<Type, Type>{
    typedef new_alloc<Type> other;
};

template <typename Type, typename Rebind>
struct new_alloc_rebind {
    typedef std::allocator<Rebind> other;
};

// designed for typed_block (to use with std::vector)
template <class Type>
class new_alloc
{
public:
    // type definitions
    typedef Type value_type;
    typedef Type* pointer;
    typedef const Type* const_pointer;
    typedef Type& reference;
    typedef const Type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    // rebind allocator to type Rebind, use new_alloc only if Rebind == Type
    template <class Rebind>
    struct rebind {
        typedef typename new_alloc_rebind<Type, Rebind>::other other;
    };

    // return address of values
    pointer address(reference value) const
    {
        return &value;
    }
    const_pointer address(const_reference value) const
    {
        return &value;
    }

    new_alloc() throw () { }
    new_alloc(const new_alloc&) throw () { }
    template <class Rebind>
    new_alloc(const new_alloc<Rebind>&) throw () { }
    ~new_alloc() throw () { }

    template <class Rebind>
    operator std::allocator<Rebind>()
    {
        static std::allocator<Rebind> helper_allocator;
        return helper_allocator;
    }

    // return maximum number of elements that can be allocated
    size_type max_size() const throw ()
    {
        return std::numeric_limits<size_type>::max() / sizeof(Type);
    }

    // allocate but don't initialize num elements of type Type
    pointer allocate(size_type num, const void* = 0)
    {
        return static_cast<Type*>(Type::operator new (num * sizeof(Type)));
    }

    // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 402. wrong new expression in [some_] allocator::construct
    // initialize elements of allocated storage p with value value
    void construct(pointer p, const Type& value)
    {
        // initialize memory with placement new
        ::new ((void*)p)Type(value);
    }

#ifdef __GXX_EXPERIMENTAL_CXX0X__
    template <typename ... Args>
    void construct(pointer p, Args&& ... args)
    {
        ::new ((void*)p)Type(std::forward<Args>(args) ...);
    }
#endif

    // destroy elements of initialized storage p
    void destroy(pointer p)
    {
        // destroy objects by calling their destructor
        p->~Type();
    }

    // deallocate storage p of deleted elements
    void deallocate(pointer p, size_type /*num*/)
    {
        Type::operator delete (p);
    }
};

// return that all specializations of this allocator are interchangeable
template <class Type1, class Type2>
inline bool operator == (const new_alloc<Type1>&,
                         const new_alloc<Type2>&) throw ()
{
    return true;
}

template <class Type1, class Type2>
inline bool operator != (const new_alloc<Type1>&,
                         const new_alloc<Type2>&) throw ()
{
    return false;
}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_NEW_ALLOC_HEADER
// vim: et:ts=4:sw=4
