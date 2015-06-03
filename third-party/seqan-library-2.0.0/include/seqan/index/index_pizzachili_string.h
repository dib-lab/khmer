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
// Author: Konrad Ludwig Moritz Rudolph <konrad.rudolph@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_STRING_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_STRING_H

#include <cstdlib> // for std::malloc, std::realloc, std::free

#include <seqan/index/pizzachili_api.h>

namespace SEQAN_NAMESPACE_MAIN {

/*!
 * @class PizzaChiliString Pizza &amp; Chili String
 *
 * @extends String
 *
 * @headerfile <seqan/index.h>
 *
 * @brief String used by the Pizza &amp; Chili indices.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class String<TValue, PizzaChili<TSpec> >;
 *
 * @tparam TSpec Tag specifying the Pizza &amp; Chili index library to use. Types:
 *               PizzaChiliIndexTags
 * @tparam TValue The value type, that is the type of them items/characters
 *                stored in the string.This type must be a simple type and it
 *                must hold that <tt>sizeof(TValue) == 1</tt>.
 *
 * The string is lazy in the sense that it holds a reference to the compressed
 * index structure it is associated with. Only when the text is actually read,
 * the index is queried for the text. If only a substring is needed, this string
 * tries to query only a substring.
 *
 * @see PizzaChiliIndex
 * @see PizzaChiliIndexTags
 */

template <typename TSpec>
struct PizzaChili;

template <typename TValue, typename TSpec>
class String<TValue, PizzaChili<TSpec> > {
public:
    typedef TValue* TIter;
    typedef TValue const* TConstIter;
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

    impl::index_t index_handle;
    bool owned;
    mutable TIter data_begin;
    mutable TIter data_end;

    String() : index_handle(0), owned(true), data_begin(0), data_end(0) { }

    template <typename TText>
    String(TText& other)
        : index_handle(0), owned(true), data_begin(0), data_end(0)
    {
SEQAN_CHECKPOINT
        assign(*this, other);
    }

    String(impl::index_t index_handle)
        : index_handle(index_handle), owned(true), data_begin(0), data_end(0) { }

    String(String& other)
        : index_handle(0), owned(true), data_begin(0), data_end(0)
    {
SEQAN_CHECKPOINT
        assign(*this, other);
    }

    String(String const& other)
        : index_handle(0), owned(true), data_begin(0), data_end(0)
    {
SEQAN_CHECKPOINT
        assign(*this, other);
    }

    String& operator =(String const& other) {
SEQAN_CHECKPOINT
        if (this != &other)
            assign(*this, other);
        return *this;
    }

    ~String() {
SEQAN_CHECKPOINT
        clear(*this);
    }

    // Convenience ...

    template <typename TPos>
    inline TValue& operator [](TPos index) {
SEQAN_CHECKPOINT
        return value(*this, index);
    }

    template <typename TPos>
    inline TValue operator [](TPos index) const {
SEQAN_CHECKPOINT
        return value(*this, index);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct IsContiguous<String<TValue, PizzaChili<TSpec> > > {
    enum { VALUE = true };
};

template <typename TValue, typename TSpec>
struct IsContiguous<String<TValue, PizzaChili<TSpec> > const> {
    enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TTag>
struct Iterator<String<TValue, PizzaChili<TSpec> >, TTag> {
    typedef typename String<TValue, PizzaChili<TSpec> >::TIter Type;
};

template <typename TValue, typename TSpec, typename TTag>
struct Iterator<String<TValue, PizzaChili<TSpec> > const, TTag> {
    typedef typename String<TValue, PizzaChili<TSpec> >::TConstIter Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Infix<String<TValue, PizzaChili<TSpec> > > {
    typedef String<TValue, PizzaChili<TSpec> > Type;
};

template <typename TValue, typename TSpec>
struct Infix<String<TValue, PizzaChili<TSpec> > const> {
    typedef String<TValue, PizzaChili<TSpec> > const Type;
};

template <typename TValue, typename TSpec>
struct Prefix<String<TValue, PizzaChili<TSpec> > >
    : Infix<String<TValue, PizzaChili<TSpec> > > { };

template <typename TValue, typename TSpec>
struct Prefix<String<TValue, PizzaChili<TSpec> > const>
    : Infix<String<TValue, PizzaChili<TSpec> > const> { };

template <typename TValue, typename TSpec>
struct Suffix<String<TValue, PizzaChili<TSpec> > >
    : Infix<String<TValue, PizzaChili<TSpec> > > { };

template <typename TValue, typename TSpec>
struct Suffix<String<TValue, PizzaChili<TSpec> > const>
    : Infix<String<TValue, PizzaChili<TSpec> > const> { };

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct DefaultOverflowImplicit<String<TValue, PizzaChili<TSpec> > > {
    typedef Exact Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(String<TValue, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    me.index_handle = 0;
    if (me.owned) {
        _deallocateStorage(me, me.data_begin, me.data_end - me.data_begin);
        me.data_begin = 0;
        me.data_end = 0;
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSource, typename TExpand>
inline void
assign(
    String<TValue, PizzaChili<TSpec> >& target,
    TSource const& source,
    Tag<TExpand>
) {
    target.owned = true;
    AssignString_<Tag<TExpand> >::assign_(target, source);
}

template<typename TValue, typename TSpec, typename TSource, typename TExpand>
inline void
assign(
    String<TValue, PizzaChili<TSpec> >& target,
    TSource const* source,
    Tag<TExpand>
) {
    target.owned = true;
    AssignString_<Tag<TExpand> >::assign_(target, source);
}

template<typename TValue, typename TSpec, typename TExpand>
inline void
assign(
    String<TValue, PizzaChili<TSpec> >& target,
    String<TValue, PizzaChili<TSpec> > const& source,
    Tag<TExpand>
) {
    target.owned = true;

    if (source.index_handle != 0) {
        clear(target);
        target.index_handle = source.index_handle;
    }
    else
        AssignString_<Tag<TExpand> >::assign_(target, source);
}

//////////////////////////////////////////////////////////////////////////////

/*
template <typename TValue, typename TSpec>
inline typename Value<String<TValue, PizzaChili<TSpec> > >::Type*
_pizzaChiliAllocate(
    String<TValue, PizzaChili<TSpec> >& me,
    typename Size<String<TValue, PizzaChili<TSpec> > >::Type new_capacity
) {
SEQAN_ CHECKPOINT
    TValue* old = me.data_begin;
    me.data_begin = static_cast<TValue*>(std::malloc(new_capacity));
    return old == me.data_begin ? 0 : old;
}
*/

template <typename TValue, typename TSpec>
inline typename Value<String<TValue, PizzaChili<TSpec> > >::Type*
_pizzaChiliReallocate(
    String<TValue, PizzaChili<TSpec> >& me,
    typename Size<String<TValue, PizzaChili<TSpec> > >::Type new_capacity
) {
SEQAN_CHECKPOINT
    if (new_capacity <= capacity(me))
        return 0;

    SEQAN_ASSERT_NEQ(me.data_begin, 0);
    //if (me.data_begin == 0)
    //    return _pizzaChiliAllocate(me, new_capacity);

    me.data_begin =
        static_cast<TValue*>(std::realloc(me.data_begin, new_capacity));
    // std::realloc does the cleanup itself.
    return 0;
}

template <typename TValue>
inline void
_pizzaChiliDeallocate(TValue* begin) {
SEQAN_CHECKPOINT
    if (begin != 0)
        std::free(begin);
}

template <typename TValue, typename TSpec>
struct AllocHelper_ {
    typedef String<TValue, PizzaChili<TSpec> > string_type;
    typedef typename Value<string_type>::Type* pointer_type;
    typedef typename Size<string_type>::Type size_type;

    /*
    static pointer_type allocate(string_type& me, size_type new_capacity) {
SEQAN_ CHECKPOINT
        return _pizzaChiliAllocate(me, new_capacity);
    }
    */

    static pointer_type reallocate(string_type& me, size_type new_capacity) {
SEQAN_CHECKPOINT
        return _pizzaChiliReallocate(me, new_capacity);
    }

    static void deallocate(pointer_type begin) {
SEQAN_CHECKPOINT
        _pizzaChiliDeallocate(begin);
    }
};

// This specialization is needed because the FM index requires than an
// overshoot be reserved after each string.
template <typename TValue>
struct AllocHelper_<TValue, PizzaChiliFM> {
    typedef PizzaChiliFM TSpec;
    typedef String<TValue, PizzaChili<TSpec> > string_type;
    typedef typename Value<string_type>::Type* pointer_type;
    typedef typename Size<string_type>::Type size_type;

    static size_type real_capacity(size_type new_capacity) {
SEQAN_CHECKPOINT
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
        // The numerical constant values are taken from the fm_build.c example
        // in th FM library folder.
        return new_capacity + TCodeProvider::init_ds_ssort(500, 2000);
    }

    /*
    static pointer_type allocate(string_type& me, size_type new_capacity) {
SEQAN_ CHECKPOINT
        return _pizzaChiliAllocate(me, real_capacity(new_capacity));
    }
    */

    static pointer_type reallocate(string_type& me, size_type new_capacity) {
SEQAN_CHECKPOINT
        return _pizzaChiliReallocate(me, real_capacity(new_capacity));
    }

    static void deallocate(pointer_type begin) {
SEQAN_CHECKPOINT
        _pizzaChiliDeallocate(begin);
    }
};

/*
template <typename TValue, typename TSpec>
inline typename Value<String<TValue, PizzaChili<TSpec> > >::Type*
_allocateStorage(
    String<TValue, PizzaChili<TSpec> >& me,
    typename Size<String<TValue, PizzaChili<TSpec> > >::Type new_capacity
) {
SEQAN_ CHECKPOINT
    return AllocHelper_<TValue, TSpec>::allocate(me, new_capacity);
}
*/

template <typename TValue, typename TSpec>
inline typename Value<String<TValue, PizzaChili<TSpec> > >::Type*
_reallocateStorage(
    String<TValue, PizzaChili<TSpec> >& me,
    typename Size<String<TValue, PizzaChili<TSpec> > >::Type new_capacity,
    Exact
) {
SEQAN_CHECKPOINT
    return AllocHelper_<TValue, TSpec>::reallocate(me, new_capacity);
}

template <typename TValue, typename TSpec>
inline void
_deallocateStorage(
   String<TValue, PizzaChili<TSpec> >& /*me*/,
   TValue* begin,
   typename Size<String<TValue, PizzaChili<TSpec> > >::Type /*count*/
) {
SEQAN_CHECKPOINT
    AllocHelper_<TValue, TSpec>::deallocate(begin);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<String<TValue, PizzaChili<TSpec> > >::Type
length(String<TValue, PizzaChili<TSpec> > const& me) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

    if (me.data_begin != 0) {
SEQAN_CHECKPOINT
        return me.data_end - me.data_begin;
    }
    else if (me.index_handle != 0) {
SEQAN_CHECKPOINT
        impl::ulong_t len;
        impl::error_t e =
            TCodeProvider::get_length(me.index_handle, &len);
        if (e != 0) {
            SEQAN_REPORT(TCodeProvider::error_index(e));
            struct { } ex;
            SEQAN_THROW(ex);
        }

        return len;
    }
    else
        return 0;
}

template <typename TValue, typename TSpec>
inline void
_setLength(
    String<TValue, PizzaChili<TSpec> >& me,
    size_t new_length
) {
SEQAN_CHECKPOINT
    me.data_end = me.data_begin + new_length;
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TValue, typename TSpec>
    inline void
    queryText(String<TValue, PizzaChili<TSpec> > const& me) {
SEQAN_CHECKPOINT
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
        if (me.data_begin != 0) {
SEQAN_CHECKPOINT
            return;
        }
        if (me.index_handle == 0) {
SEQAN_CHECKPOINT
            me.data_begin = static_cast<TValue*>(std::malloc(1));
            me.data_begin[0] = '\0';
            me.data_end = me.data_begin;
        }
        else {
SEQAN_CHECKPOINT
            impl::uchar_t* snippet;
            impl::ulong_t len;
            impl::error_t e =
                TCodeProvider::extract(
                    me.index_handle,
                    0,
                    length(me) - 1,
                    &snippet,
                    &len
                );

            if (e != 0) {
                SEQAN_REPORT(TCodeProvider::error_index(e));
                struct { } ex;
                SEQAN_THROW(ex);
            }

            me.data_begin = reinterpret_cast<TValue*>(snippet);
            me.data_end = me.data_begin + len;
        }
    }
} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> >, Tag<TSpec> const>::Type
begin(
    String<TValue, PizzaChili<TSpec> >& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_begin;
}

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> > const, Tag<TSpec> const>::Type
begin(
    String<TValue, PizzaChili<TSpec> > const& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_begin;
}

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> >, Tag<TSpec> const>::Type
end(
    String<TValue, PizzaChili<TSpec> >& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_end;
}

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> > const, Tag<TSpec> const>::Type
end(
    String<TValue, PizzaChili<TSpec> > const& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_end;
}

//////////////////////////////////////////////////////////////////////////////

// NOTE: The following might seem redundant -- and it is -- but I don't see an
// easier way to implement it. The helper constructs the infix, prefix and
// suffix for a given Pizza & Chili string. The partial specializations cover
// the following cases:
//
//  - integral positions are given
//  - positions are iterators
//  - positions are const iterators
//
// Additionally, every of those cases handles const and non-const string
// arguments, for all three algorithms. Thus, there's a total of 3 * 3 * 2
// cases, resulting in 18 functions.

namespace impl {
    template <typename TValue, typename TSpec, typename TPos>
    struct substringHelperPizzaChili {
        typedef String<TValue, PizzaChili<TSpec> > TString;
        typedef typename Infix<TString>::Type TResult;
        typedef typename Infix<TString const>::Type TConstResult;
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

        static inline TResult
        infix(TString& me, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            TResult ret;

            if (me.data_begin != 0) {
SEQAN_CHECKPOINT
                ret.owned = false;
                ret.data_begin = me.data_begin + begin;
                ret.data_end = me.data_begin + end;
            }
            else if (me.index_handle != 0) {
SEQAN_CHECKPOINT
                impl::uchar_t* snippet;
                impl::ulong_t len;
                impl::error_t e =
                    TCodeProvider::extract(
                        me.index_handle,
                        begin,
                        end - 1,
                        &snippet,
                        &len
                    );

                if (e != 0) {
                    SEQAN_REPORT(TCodeProvider::error_index(e));
                    struct { } ex;
                    SEQAN_THROW(ex);
                }

                ret.data_begin = reinterpret_cast<TValue*>(snippet);
                ret.data_end = ret.data_begin + len;
            }

            return ret;
        }

        static inline TConstResult
        infix(TString const& me, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            // This cast is safe as the content is never written to.
            return infix(const_cast<TString&>(me), begin, end);
        }

        static inline TResult
        prefix(TString& me, TPos end) {
SEQAN_CHECKPOINT
            return infix(me, 0, end);
        }

        static inline TConstResult
        prefix(TString const& me, TPos end) {
SEQAN_CHECKPOINT
            return infix(me, 0, end);
        }

        static inline TResult
        suffix(TString& me, TPos begin) {
SEQAN_CHECKPOINT
            return infix(me, begin, length(me));
        }

        static inline TConstResult
        suffix(TString const& me, TPos begin) {
SEQAN_CHECKPOINT
            return infix(me, begin, length(me));
        }
    }; // substringHelperPizzaChili

    template <typename TValue, typename TSpec>
    struct substringHelperPizzaChili<
        TValue,
        TSpec,
        typename Iterator<
            String<TValue, PizzaChili<TSpec> >,
            typename DefaultIteratorSpec<String<TValue, PizzaChili<TSpec> > >::Type
        >::Type
    > {
        typedef String<TValue, PizzaChili<TSpec> > TString;
        typedef
            typename Iterator<
                TString, typename DefaultIteratorSpec<TString>::Type
            >::Type TPos;
        typedef typename Infix<TString>::Type TResult;
        typedef typename Infix<TString const>::Type TConstResult;

        static inline TResult
        infix(TString& /*me*/, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            // Iterators were used, therefore it's safe to assume that the
            // string is already in memory.
            TResult ret;
            ret.owned = false;
            ret.data_begin = begin;
            ret.data_end = end;
            return ret;
        }

        static inline TConstResult
        infix(TString const& me, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            // This cast is safe as the content is never written to.
            return infix(const_cast<TString&>(me), begin, end);
        }

        static inline TResult
        prefix(TString& me, TPos end) {
SEQAN_CHECKPOINT
            return infix(me, begin(me), end);
        }

        static inline TConstResult
        prefix(TString const& me, TPos end) {
SEQAN_CHECKPOINT
            // This cast is safe as the content is never written to.
            return infix(me, begin(const_cast<TString&>(me)), end);
        }

        static inline TResult
        suffix(TString& me, TPos begin) {
SEQAN_CHECKPOINT
            return infix(me, begin, end(me));
        }

        static inline TConstResult
        suffix(TString const& me, TPos begin) {
SEQAN_CHECKPOINT
            // This cast is safe as the content is never written to.
            return infix(me, begin, end(const_cast<TString&>(me)));
        }
    }; // substringHelperPizzaChili

    // Specialization for const iterators.
    template <typename TValue, typename TSpec>
    struct substringHelperPizzaChili<
        TValue,
        TSpec,
        typename Iterator<
            String<TValue, PizzaChili<TSpec> > const,
            typename DefaultIteratorSpec<String<TValue, PizzaChili<TSpec> > >::Type
        >::Type
    > : substringHelperPizzaChili<
        TValue, TSpec,
        typename Iterator<
            String<TValue, PizzaChili<TSpec> >,
            typename DefaultIteratorSpec<String<TValue, PizzaChili<TSpec> > >::Type
        >::Type
    > {
        typedef String<TValue, PizzaChili<TSpec> > TString;
        typedef
            typename Iterator<
                TString const,
                typename DefaultIteratorSpec<TString>::Type
            >::Type TPos;
        typedef typename Infix<TString const>::Type TResult;
        typedef typename Infix<TString const>::Type TConstResult;
        typedef substringHelperPizzaChili<
            TValue, TSpec,
            typename Iterator<
                TString,
                typename DefaultIteratorSpec<TString>::Type
            >::Type
        > TBase;

        static inline TConstResult
        infix(TString const& me, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            // This cast is safe as the content is never written to.
            return TBase::infix(
                me,
                const_cast<typename TBase::TPos>(begin),
                const_cast<typename TBase::TPos>(end)
            );
        }

        static inline TConstResult
        prefix(TString const& me, TPos end) {
SEQAN_CHECKPOINT
            return infix(me, begin(me), end);
        }

        static inline TConstResult
        suffix(TString const& me, TPos begin) {
SEQAN_CHECKPOINT
            return infix(me, begin, end(me));
        }
    };

} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<String<TValue, PizzaChili<TSpec> > >::Type
infix(
    String<TValue, PizzaChili<TSpec> > const& me,
    TPosBegin begin,
    TPosEnd end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPosBegin>::infix(me, begin, end);
}

template <typename TValue, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<String<TValue, PizzaChili<TSpec> > >::Type
infix(
    String<TValue, PizzaChili<TSpec> >& me,
    TPosBegin begin,
    TPosEnd end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPosBegin>::infix(me, begin, end);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline typename Prefix<String<TValue, PizzaChili<TSpec> > >::Type
prefix(
    String<TValue, PizzaChili<TSpec> > const& me,
    TPos end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::prefix(me, end);
}

template <typename TValue, typename TSpec, typename TPos>
inline typename Prefix<String<TValue, PizzaChili<TSpec> > >::Type
prefix(
    String<TValue, PizzaChili<TSpec> >& me,
    TPos end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::prefix(me, end);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline typename Suffix<String<TValue, PizzaChili<TSpec> > >::Type
suffix(
    String<TValue, PizzaChili<TSpec> > const& me,
    TPos begin
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::suffix(me, begin);
}

template <typename TValue, typename TSpec, typename TPos>
inline typename Suffix<String<TValue, PizzaChili<TSpec> > >::Type
suffix(
    String<TValue, PizzaChili<TSpec> >& me,
    TPos begin
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::suffix(me, begin);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_STRING_H
