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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Superclasses to extend a class by standard operators that make use of a
// minimal set of basic operators.
// Taken from the Math Library in Boost version 1.47.
// ==========================================================================

#ifndef SEQAN_MATH_OPERATORS_H_
#define SEQAN_MATH_OPERATORS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

namespace detail
{

    template <typename T> class empty_base {};

} // namespace detail

// In this section we supply the xxxx1 and xxxx2 forms of the operator
// templates, which are explicitly targeted at the 1-type-argument and
// 2-type-argument operator forms, respectively. Some compilers get confused
// when inline friend functions are overloaded in namespaces other than the
// global namespace. When SEQAN_NO_OPERATORS_IN_NAMESPACE is defined, all of
// these templates must go in the global namespace.

//  Basic operator classes (contributed by Dave Abrahams) ------------------//

//  Note that friend functions defined in a class are implicitly inline.
//  See the C++ std, 11.4 [class.friend] paragraph 5

template <class T, class U, class B = detail::empty_base<T> >
struct less_than_comparable2 : B
{
     friend bool operator<=(const T& x, const U& y) { return !static_cast<bool>(x > y); }
     friend bool operator>=(const T& x, const U& y) { return !static_cast<bool>(x < y); }
     friend bool operator>(const U& x, const T& y)  { return y < x; }
     friend bool operator<(const U& x, const T& y)  { return y > x; }
     friend bool operator<=(const U& x, const T& y) { return !static_cast<bool>(y < x); }
     friend bool operator>=(const U& x, const T& y) { return !static_cast<bool>(y > x); }
};

template <class T, class B = detail::empty_base<T> >
struct less_than_comparable1 : B
{
     friend bool operator>(const T& x, const T& y)  { return y < x; }
     friend bool operator<=(const T& x, const T& y) { return !static_cast<bool>(y < x); }
     friend bool operator>=(const T& x, const T& y) { return !static_cast<bool>(x < y); }
};

template <class T, class U, class B = detail::empty_base<T> >
struct equality_comparable2 : B
{
     friend bool operator==(const U& y, const T& x) { return x == y; }
     friend bool operator!=(const U& y, const T& x) { return !static_cast<bool>(x == y); }
     friend bool operator!=(const T& y, const U& x) { return !static_cast<bool>(y == x); }
};

template <class T, class B = detail::empty_base<T> >
struct equality_comparable1 : B
{
     friend bool operator!=(const T& x, const T& y) { return !static_cast<bool>(x == y); }
};

// A macro which produces "name_2left" from "name".
#define SEQAN_OPERATOR2_LEFT(name) name##2##_##left

//  NRVO-friendly implementation (contributed by Daniel Frey) ---------------//

#if defined(SEQAN_HAS_NRVO) || defined(SEQAN_FORCE_SYMMETRIC_OPERATORS)

// This is the optimal implementation for ISO/ANSI C++,
// but it requires the compiler to implement the NRVO.
// If the compiler has no NRVO, this is the best symmetric
// implementation available.

#define SEQAN_BINARY_OPERATOR_COMMUTATIVE( NAME, OP )                         \
template <class T, class U, class B = detail::empty_base<T> >        \
struct NAME##2 : B                                                            \
{                                                                             \
  friend T operator OP( const T& lhs, const U& rhs )                          \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                              \
  friend T operator OP( const U& lhs, const T& rhs )                          \
    { T nrv( rhs ); nrv OP##= lhs; return nrv; }                              \
};                                                                            \
                                                                              \
template <class T, class B = detail::empty_base<T> >                 \
struct NAME##1 : B                                                            \
{                                                                             \
  friend T operator OP( const T& lhs, const T& rhs )                          \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                              \
};

#define SEQAN_BINARY_OPERATOR_NON_COMMUTATIVE( NAME, OP )               \
template <class T, class U, class B = detail::empty_base<T> >  \
struct NAME##2 : B                                                      \
{                                                                       \
  friend T operator OP( const T& lhs, const U& rhs )                    \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                        \
};                                                                      \
                                                                        \
template <class T, class U, class B = detail::empty_base<T> >  \
struct SEQAN_OPERATOR2_LEFT(NAME) : B                                   \
{                                                                       \
  friend T operator OP( const U& lhs, const T& rhs )                    \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                        \
};                                                                      \
                                                                        \
template <class T, class B = detail::empty_base<T> >           \
struct NAME##1 : B                                                      \
{                                                                       \
  friend T operator OP( const T& lhs, const T& rhs )                    \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                        \
};

#else // defined(SEQAN_HAS_NRVO) || defined(SEQAN_FORCE_SYMMETRIC_OPERATORS)

// For compilers without NRVO the following code is optimal, but not
// symmetric!  Note that the implementation of
// SEQAN_OPERATOR2_LEFT(NAME) only looks cool, but doesn't provide
// optimization opportunities to the compiler :)

#define SEQAN_BINARY_OPERATOR_COMMUTATIVE( NAME, OP )                   \
template <class T, class U, class B = detail::empty_base<T> >  \
struct NAME##2 : B                                                      \
{                                                                       \
  friend T operator OP( T lhs, const U& rhs ) { return lhs OP##= rhs; } \
  friend T operator OP( const U& lhs, T rhs ) { return rhs OP##= lhs; } \
};                                                                      \
                                                                        \
template <class T, class B = detail::empty_base<T> >           \
struct NAME##1 : B                                                      \
{                                                                       \
  friend T operator OP( T lhs, const T& rhs ) { return lhs OP##= rhs; } \
};

#define SEQAN_BINARY_OPERATOR_NON_COMMUTATIVE( NAME, OP )               \
template <class T, class U, class B = detail::empty_base<T> >  \
struct NAME##2 : B                                                      \
{                                                                       \
  friend T operator OP( T lhs, const U& rhs ) { return lhs OP##= rhs; } \
};                                                                      \
                                                                        \
template <class T, class U, class B = detail::empty_base<T> >  \
struct SEQAN_OPERATOR2_LEFT(NAME) : B                                   \
{                                                                       \
  friend T operator OP( const U& lhs, const T& rhs )                    \
    { return T( lhs ) OP##= rhs; }                                      \
};                                                                      \
                                                                        \
template <class T, class B = detail::empty_base<T> >           \
struct NAME##1 : B                                                      \
{                                                                       \
  friend T operator OP( T lhs, const T& rhs ) { return lhs OP##= rhs; } \
};

#endif // defined(SEQAN_HAS_NRVO) || defined(SEQAN_FORCE_SYMMETRIC_OPERATORS)

SEQAN_BINARY_OPERATOR_COMMUTATIVE( multipliable, * )
SEQAN_BINARY_OPERATOR_COMMUTATIVE( addable, + )
SEQAN_BINARY_OPERATOR_NON_COMMUTATIVE( subtractable, - )
SEQAN_BINARY_OPERATOR_NON_COMMUTATIVE( dividable, / )
SEQAN_BINARY_OPERATOR_NON_COMMUTATIVE( modable, % )
SEQAN_BINARY_OPERATOR_COMMUTATIVE( xorable, ^ )
SEQAN_BINARY_OPERATOR_COMMUTATIVE( andable, & )
SEQAN_BINARY_OPERATOR_COMMUTATIVE( orable, | )

#undef SEQAN_BINARY_OPERATOR_COMMUTATIVE
#undef SEQAN_BINARY_OPERATOR_NON_COMMUTATIVE
#undef SEQAN_OPERATOR2_LEFT

//  incrementable and decrementable contributed by Jeremy Siek

template <class T, class B = detail::empty_base<T> >
struct incrementable : B
{
  friend T operator++(T& x, int)
  {
    incrementable_type nrv(x);
    ++x;
    return nrv;
  }
private: // The use of this typedef works around a Borland bug
  typedef T incrementable_type;
};

template <class T, class B = detail::empty_base<T> >
struct decrementable : B
{
  friend T operator--(T& x, int)
  {
    decrementable_type nrv(x);
    --x;
    return nrv;
  }
private: // The use of this typedef works around a Borland bug
  typedef T decrementable_type;
};

//  Iterator operator classes (contributed by Jeremy Siek) ------------------//

template <class T, class P, class B = detail::empty_base<T> >
struct dereferenceable : B
{
  P operator->() const
  {
    return &*static_cast<const T&>(*this);
  }
};

template <class T, class I, class R, class B = detail::empty_base<T> >
struct indexable : B
{
  R operator[](I n) const
  {
    return *(static_cast<const T&>(*this) + n);
  }
};

//  More operator classes (contributed by Daryle Walker) --------------------//
//  (NRVO-friendly implementation contributed by Daniel Frey) ---------------//

#if defined(SEQAN_HAS_NRVO) || defined(SEQAN_FORCE_SYMMETRIC_OPERATORS)

#define SEQAN_BINARY_OPERATOR( NAME, OP )                                     \
template <class T, class U, class B = detail::empty_base<T> >        \
struct NAME##2 : B                                                            \
{                                                                             \
  friend T operator OP( const T& lhs, const U& rhs )                          \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                              \
};                                                                            \
                                                                              \
template <class T, class B = detail::empty_base<T> >                 \
struct NAME##1 : B                                                            \
{                                                                             \
  friend T operator OP( const T& lhs, const T& rhs )                          \
    { T nrv( lhs ); nrv OP##= rhs; return nrv; }                              \
};

#else // defined(SEQAN_HAS_NRVO) || defined(SEQAN_FORCE_SYMMETRIC_OPERATORS)

#define SEQAN_BINARY_OPERATOR( NAME, OP )                                     \
template <class T, class U, class B = detail::empty_base<T> >        \
struct NAME##2 : B                                                            \
{                                                                             \
  friend T operator OP( T lhs, const U& rhs ) { return lhs OP##= rhs; }       \
};                                                                            \
                                                                              \
template <class T, class B = detail::empty_base<T> >                 \
struct NAME##1 : B                                                            \
{                                                                             \
  friend T operator OP( T lhs, const T& rhs ) { return lhs OP##= rhs; }       \
};

#endif // defined(SEQAN_HAS_NRVO) || defined(SEQAN_FORCE_SYMMETRIC_OPERATORS)

SEQAN_BINARY_OPERATOR( left_shiftable, << )
SEQAN_BINARY_OPERATOR( right_shiftable, >> )

#undef SEQAN_BINARY_OPERATOR

template <class T, class U, class B = detail::empty_base<T> >
struct equivalent2 : B
{
  friend bool operator==(const T& x, const U& y)
  {
    return !static_cast<bool>(x < y) && !static_cast<bool>(x > y);
  }
};

template <class T, class B = detail::empty_base<T> >
struct equivalent1 : B
{
  friend bool operator==(const T&x, const T&y)
  {
    return !static_cast<bool>(x < y) && !static_cast<bool>(y < x);
  }
};

template <class T, class U, class B = detail::empty_base<T> >
struct partially_ordered2 : B
{
  friend bool operator<=(const T& x, const U& y)
    { return static_cast<bool>(x < y) || static_cast<bool>(x == y); }
  friend bool operator>=(const T& x, const U& y)
    { return static_cast<bool>(x > y) || static_cast<bool>(x == y); }
  friend bool operator>(const U& x, const T& y)
    { return y < x; }
  friend bool operator<(const U& x, const T& y)
    { return y > x; }
  friend bool operator<=(const U& x, const T& y)
    { return static_cast<bool>(y > x) || static_cast<bool>(y == x); }
  friend bool operator>=(const U& x, const T& y)
    { return static_cast<bool>(y < x) || static_cast<bool>(y == x); }
};

template <class T, class B = detail::empty_base<T> >
struct partially_ordered1 : B
{
  friend bool operator>(const T& x, const T& y)
    { return y < x; }
  friend bool operator<=(const T& x, const T& y)
    { return static_cast<bool>(x < y) || static_cast<bool>(x == y); }
  friend bool operator>=(const T& x, const T& y)
    { return static_cast<bool>(y < x) || static_cast<bool>(x == y); }
};

//  Combined operator classes (contributed by Daryle Walker) ----------------//

template <class T, class U, class B = detail::empty_base<T> >
struct totally_ordered2
    : less_than_comparable2<T, U
    , equality_comparable2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct totally_ordered1
    : less_than_comparable1<T
    , equality_comparable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct additive2
    : addable2<T, U
    , subtractable2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct additive1
    : addable1<T
    , subtractable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct multiplicative2
    : multipliable2<T, U
    , dividable2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct multiplicative1
    : multipliable1<T
    , dividable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct integer_multiplicative2
    : multiplicative2<T, U
    , modable2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct integer_multiplicative1
    : multiplicative1<T
    , modable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct arithmetic2
    : additive2<T, U
    , multiplicative2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct arithmetic1
    : additive1<T
    , multiplicative1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct integer_arithmetic2
    : additive2<T, U
    , integer_multiplicative2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct integer_arithmetic1
    : additive1<T
    , integer_multiplicative1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct bitwise2
    : xorable2<T, U
    , andable2<T, U
    , orable2<T, U, B
      > > > {};

template <class T, class B = detail::empty_base<T> >
struct bitwise1
    : xorable1<T
    , andable1<T
    , orable1<T, B
      > > > {};

template <class T, class B = detail::empty_base<T> >
struct unit_steppable
    : incrementable<T
    , decrementable<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct shiftable2
    : left_shiftable2<T, U
    , right_shiftable2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct shiftable1
    : left_shiftable1<T
    , right_shiftable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct ring_operators2
    : additive2<T, U
    , subtractable2_left<T, U
    , multipliable2<T, U, B
      > > > {};

template <class T, class B = detail::empty_base<T> >
struct ring_operators1
    : additive1<T
    , multipliable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct ordered_ring_operators2
    : ring_operators2<T, U
    , totally_ordered2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct ordered_ring_operators1
    : ring_operators1<T
    , totally_ordered1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct field_operators2
    : ring_operators2<T, U
    , dividable2<T, U
    , dividable2_left<T, U, B
      > > > {};

template <class T, class B = detail::empty_base<T> >
struct field_operators1
    : ring_operators1<T
    , dividable1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct ordered_field_operators2
    : field_operators2<T, U
    , totally_ordered2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct ordered_field_operators1
    : field_operators1<T
    , totally_ordered1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct euclidian_ring_operators2
    : ring_operators2<T, U
    , dividable2<T, U
    , dividable2_left<T, U
    , modable2<T, U
    , modable2_left<T, U, B
      > > > > > {};

template <class T, class B = detail::empty_base<T> >
struct euclidian_ring_operators1
    : ring_operators1<T
    , dividable1<T
    , modable1<T, B
      > > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct ordered_euclidian_ring_operators2
    : totally_ordered2<T, U
    , euclidian_ring_operators2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct ordered_euclidian_ring_operators1
    : totally_ordered1<T
    , euclidian_ring_operators1<T, B
      > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct euclidean_ring_operators2
    : ring_operators2<T, U
    , dividable2<T, U
    , dividable2_left<T, U
    , modable2<T, U
    , modable2_left<T, U, B
      > > > > > {};

template <class T, class B = detail::empty_base<T> >
struct euclidean_ring_operators1
    : ring_operators1<T
    , dividable1<T
    , modable1<T, B
      > > > {};

template <class T, class U, class B = detail::empty_base<T> >
struct ordered_euclidean_ring_operators2
    : totally_ordered2<T, U
    , euclidean_ring_operators2<T, U, B
      > > {};

template <class T, class B = detail::empty_base<T> >
struct ordered_euclidean_ring_operators1
    : totally_ordered1<T
    , euclidean_ring_operators1<T, B
      > > {};

template <class T, class P, class B = detail::empty_base<T> >
struct input_iteratable
    : equality_comparable1<T
    , incrementable<T
    , dereferenceable<T, P, B
      > > > {};

template <class T, class B = detail::empty_base<T> >
struct output_iteratable
    : incrementable<T, B
      > {};

template <class T, class P, class B = detail::empty_base<T> >
struct forward_iteratable
    : input_iteratable<T, P, B
      > {};

template <class T, class P, class B = detail::empty_base<T> >
struct bidirectional_iteratable
    : forward_iteratable<T, P
    , decrementable<T, B
      > > {};

//  To avoid repeated derivation from equality_comparable,
//  which is an indirect base class of bidirectional_iterable,
//  random_access_iteratable must not be derived from totally_ordered1
//  but from less_than_comparable1 only. (Helmut Zeisel, 02-Dec-2001)
template <class T, class P, class D, class R, class B = detail::empty_base<T> >
struct random_access_iteratable
    : bidirectional_iteratable<T, P
    , less_than_comparable1<T
    , additive2<T, D
    , indexable<T, D, R, B
      > > > > {};


// SEQAN_IMPORT_TEMPLATE1 .. SEQAN_IMPORT_TEMPLATE4 -
//
// When SEQAN_NO_OPERATORS_IN_NAMESPACE is defined we need a way to import an
// operator template into the boost namespace. SEQAN_IMPORT_TEMPLATE1 is used
// for one-argument forms of operator templates; SEQAN_IMPORT_TEMPLATE2 for
// two-argument forms. Note that these macros expect to be invoked from within
// boost.

#ifndef SEQAN_NO_OPERATORS_IN_NAMESPACE

  // The template is already in boost so we have nothing to do.
# define SEQAN_IMPORT_TEMPLATE4(template_name)
# define SEQAN_IMPORT_TEMPLATE3(template_name)
# define SEQAN_IMPORT_TEMPLATE2(template_name)
# define SEQAN_IMPORT_TEMPLATE1(template_name)

#else // SEQAN_NO_OPERATORS_IN_NAMESPACE

#  ifndef SEQAN_NO_USING_TEMPLATE

     // Bring the names in with a using-declaration
     // to avoid stressing the compiler.
#    define SEQAN_IMPORT_TEMPLATE4(template_name) using ::template_name;
#    define SEQAN_IMPORT_TEMPLATE3(template_name) using ::template_name;
#    define SEQAN_IMPORT_TEMPLATE2(template_name) using ::template_name;
#    define SEQAN_IMPORT_TEMPLATE1(template_name) using ::template_name;

#  else

     // Otherwise, because a Borland C++ 5.5 bug prevents a using declaration
     // from working, we are forced to use inheritance for that compiler.
#    define SEQAN_IMPORT_TEMPLATE4(template_name)                                             \
     template <class T, class U, class V, class W, class B = detail::empty_base<T> > \
     struct template_name : ::template_name<T, U, V, W, B> {};

#    define SEQAN_IMPORT_TEMPLATE3(template_name)                                    \
     template <class T, class U, class V, class B = detail::empty_base<T> > \
     struct template_name : ::template_name<T, U, V, B> {};

#    define SEQAN_IMPORT_TEMPLATE2(template_name)                           \
     template <class T, class U, class B = detail::empty_base<T> > \
     struct template_name : ::template_name<T, U, B> {};

#    define SEQAN_IMPORT_TEMPLATE1(template_name)                  \
     template <class T, class B = detail::empty_base<T> > \
     struct template_name : ::template_name<T, B> {};

#  endif // SEQAN_NO_USING_TEMPLATE

#endif // SEQAN_NO_OPERATORS_IN_NAMESPACE

//
// Here's where we put it all together, defining the xxxx forms of the templates
// in namespace boost. We also define specializations of is_chained_base<> for
// the xxxx, xxxx1, and xxxx2 templates, importing them into boost:: as
// necessary.
//
#ifndef SEQAN_NO_TEMPLATE_PARTIAL_SPECIALIZATION

// is_chained_base<> - a traits class used to distinguish whether an operator
// template argument is being used for base class chaining, or is specifying a
// 2nd argument type.

// A type parameter is used instead of a plain bool because Borland's compiler
// didn't cope well with the more obvious non-type template parameter.
namespace detail {
  struct true_t {};
  struct false_t {};
} // namespace detail

// Unspecialized version assumes that most types are not being used for base
// class chaining. We specialize for the operator templates defined in this
// library.
template<class T> struct is_chained_base {
  typedef detail::false_t value;
};

// Import a 4-type-argument operator template into boost (if necessary) and
// provide a specialization of 'is_chained_base<>' for it.
# define SEQAN_OPERATOR_TEMPLATE4(template_name4)                     \
  SEQAN_IMPORT_TEMPLATE4(template_name4)                              \
  template<class T, class U, class V, class W, class B>               \
  struct is_chained_base< template_name4<T, U, V, W, B> > {  \
    typedef detail::true_t value;                            \
  };

// Import a 3-type-argument operator template into boost (if necessary) and
// provide a specialization of 'is_chained_base<>' for it.
# define SEQAN_OPERATOR_TEMPLATE3(template_name3)                     \
  SEQAN_IMPORT_TEMPLATE3(template_name3)                              \
  template<class T, class U, class V, class B>                        \
  struct is_chained_base< template_name3<T, U, V, B> > {     \
    typedef detail::true_t value;                            \
  };

// Import a 2-type-argument operator template into boost (if necessary) and
// provide a specialization of 'is_chained_base<>' for it.
# define SEQAN_OPERATOR_TEMPLATE2(template_name2)                  \
  SEQAN_IMPORT_TEMPLATE2(template_name2)                           \
  template<class T, class U, class B>                              \
  struct is_chained_base< template_name2<T, U, B> > {     \
    typedef detail::true_t value;                         \
  };

// Import a 1-type-argument operator template into boost (if necessary) and
// provide a specialization of 'is_chained_base<>' for it.
# define SEQAN_OPERATOR_TEMPLATE1(template_name1)                  \
  SEQAN_IMPORT_TEMPLATE1(template_name1)                           \
  template<class T, class B>                                       \
  struct is_chained_base< template_name1<T, B> > {        \
    typedef detail::true_t value;                         \
  };

// SEQAN_OPERATOR_TEMPLATE(template_name) defines template_name<> such that it
// can be used for specifying both 1-argument and 2-argument forms. Requires the
// existence of two previously defined class templates named '<template_name>1'
// and '<template_name>2' which must implement the corresponding 1- and 2-
// argument forms.
//
// The template type parameter O == is_chained_base<U>::value is used to
// distinguish whether the 2nd argument to <template_name> is being used for
// base class chaining from another boost operator template or is describing a
// 2nd operand type. O == true_t only when U is actually an another operator
// template from the library. Partial specialization is used to select an
// implementation in terms of either '<template_name>1' or '<template_name>2'.
//

# define SEQAN_OPERATOR_TEMPLATE(template_name)                    \
template <class T                                                  \
         ,class U = T                                              \
         ,class B = detail::empty_base<T>                 \
         ,class O = typename is_chained_base<U>::value             \
         >                                                         \
struct template_name : template_name##2<T, U, B> {};               \
                                                                   \
template<class T, class U, class B>                                \
struct template_name<T, U, B, detail::true_t>             \
  : template_name##1<T, U> {};                                     \
                                                                   \
template <class T, class B>                                        \
struct template_name<T, T, B, detail::false_t>            \
  : template_name##1<T, B> {};                                     \
                                                                   \
template<class T, class U, class B, class O>                       \
struct is_chained_base< template_name<T, U, B, O> > {     \
  typedef detail::true_t value;                           \
};                                                                 \
                                                                   \
SEQAN_OPERATOR_TEMPLATE2(template_name##2)                         \
SEQAN_OPERATOR_TEMPLATE1(template_name##1)


#else // SEQAN_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#  define SEQAN_OPERATOR_TEMPLATE4(template_name4) \
        SEQAN_IMPORT_TEMPLATE4(template_name4)
#  define SEQAN_OPERATOR_TEMPLATE3(template_name3) \
        SEQAN_IMPORT_TEMPLATE3(template_name3)
#  define SEQAN_OPERATOR_TEMPLATE2(template_name2) \
        SEQAN_IMPORT_TEMPLATE2(template_name2)
#  define SEQAN_OPERATOR_TEMPLATE1(template_name1) \
        SEQAN_IMPORT_TEMPLATE1(template_name1)

   // In this case we can only assume that template_name<> is equivalent to the
   // more commonly needed template_name1<> form.
#  define SEQAN_OPERATOR_TEMPLATE(template_name)                   \
   template <class T, class B = detail::empty_base<T> >   \
   struct template_name : template_name##1<T, B> {};

#endif // SEQAN_NO_TEMPLATE_PARTIAL_SPECIALIZATION

SEQAN_OPERATOR_TEMPLATE(less_than_comparable)
SEQAN_OPERATOR_TEMPLATE(equality_comparable)
SEQAN_OPERATOR_TEMPLATE(multipliable)
SEQAN_OPERATOR_TEMPLATE(addable)
SEQAN_OPERATOR_TEMPLATE(subtractable)
SEQAN_OPERATOR_TEMPLATE2(subtractable2_left)
SEQAN_OPERATOR_TEMPLATE(dividable)
SEQAN_OPERATOR_TEMPLATE2(dividable2_left)
SEQAN_OPERATOR_TEMPLATE(modable)
SEQAN_OPERATOR_TEMPLATE2(modable2_left)
SEQAN_OPERATOR_TEMPLATE(xorable)
SEQAN_OPERATOR_TEMPLATE(andable)
SEQAN_OPERATOR_TEMPLATE(orable)

SEQAN_OPERATOR_TEMPLATE1(incrementable)
SEQAN_OPERATOR_TEMPLATE1(decrementable)

SEQAN_OPERATOR_TEMPLATE2(dereferenceable)
SEQAN_OPERATOR_TEMPLATE3(indexable)

SEQAN_OPERATOR_TEMPLATE(left_shiftable)
SEQAN_OPERATOR_TEMPLATE(right_shiftable)
SEQAN_OPERATOR_TEMPLATE(equivalent)
SEQAN_OPERATOR_TEMPLATE(partially_ordered)

SEQAN_OPERATOR_TEMPLATE(totally_ordered)
SEQAN_OPERATOR_TEMPLATE(additive)
SEQAN_OPERATOR_TEMPLATE(multiplicative)
SEQAN_OPERATOR_TEMPLATE(integer_multiplicative)
SEQAN_OPERATOR_TEMPLATE(arithmetic)
SEQAN_OPERATOR_TEMPLATE(integer_arithmetic)
SEQAN_OPERATOR_TEMPLATE(bitwise)
SEQAN_OPERATOR_TEMPLATE1(unit_steppable)
SEQAN_OPERATOR_TEMPLATE(shiftable)
SEQAN_OPERATOR_TEMPLATE(ring_operators)
SEQAN_OPERATOR_TEMPLATE(ordered_ring_operators)
SEQAN_OPERATOR_TEMPLATE(field_operators)
SEQAN_OPERATOR_TEMPLATE(ordered_field_operators)
SEQAN_OPERATOR_TEMPLATE(euclidian_ring_operators)
SEQAN_OPERATOR_TEMPLATE(ordered_euclidian_ring_operators)
SEQAN_OPERATOR_TEMPLATE(euclidean_ring_operators)
SEQAN_OPERATOR_TEMPLATE(ordered_euclidean_ring_operators)
SEQAN_OPERATOR_TEMPLATE2(input_iteratable)
SEQAN_OPERATOR_TEMPLATE1(output_iteratable)
SEQAN_OPERATOR_TEMPLATE2(forward_iteratable)
SEQAN_OPERATOR_TEMPLATE2(bidirectional_iteratable)
SEQAN_OPERATOR_TEMPLATE4(random_access_iteratable)

#undef SEQAN_OPERATOR_TEMPLATE
#undef SEQAN_OPERATOR_TEMPLATE4
#undef SEQAN_OPERATOR_TEMPLATE3
#undef SEQAN_OPERATOR_TEMPLATE2
#undef SEQAN_OPERATOR_TEMPLATE1
#undef SEQAN_IMPORT_TEMPLATE1
#undef SEQAN_IMPORT_TEMPLATE2
#undef SEQAN_IMPORT_TEMPLATE3
#undef SEQAN_IMPORT_TEMPLATE4

// The following 'operators' classes can only be used portably if the derived class
// declares ALL of the required member operators.
template <class T, class U>
struct operators2
    : totally_ordered2<T,U
    , integer_arithmetic2<T,U
    , bitwise2<T,U
      > > > {};

#ifndef SEQAN_NO_TEMPLATE_PARTIAL_SPECIALIZATION
template <class T, class U = T>
struct operators : operators2<T, U> {};

template <class T> struct operators<T, T>
#else
template <class T> struct operators
#endif
    : totally_ordered<T
    , integer_arithmetic<T
    , bitwise<T
    , unit_steppable<T
      > > > > {};

/*

//  Iterator helper classes (contributed by Jeremy Siek) -------------------//
//  (Input and output iterator helpers contributed by Daryle Walker) -------//
//  (Changed to use combined operator classes by Daryle Walker) ------------//
template <class T,
          class V,
          class D = std::ptrdiff_t,
          class P = V const *,
          class R = V const &>
struct input_iterator_helper
  : input_iteratable<T, P
  , boost::iterator<std::input_iterator_tag, V, D, P, R
    > > {};

template<class T>
struct output_iterator_helper
  : output_iteratable<T
  , boost::iterator<std::output_iterator_tag, void, void, void, void
  > >
{
  T& operator*()  { return static_cast<T&>(*this); }
  T& operator++() { return static_cast<T&>(*this); }
};

template <class T,
          class V,
          class D = std::ptrdiff_t,
          class P = V*,
          class R = V&>
struct forward_iterator_helper
  : forward_iteratable<T, P
  , boost::iterator<std::forward_iterator_tag, V, D, P, R
    > > {};

template <class T,
          class V,
          class D = std::ptrdiff_t,
          class P = V*,
          class R = V&>
struct bidirectional_iterator_helper
  : bidirectional_iteratable<T, P
  , boost::iterator<std::bidirectional_iterator_tag, V, D, P, R
    > > {};

template <class T,
          class V,
          class D = std::ptrdiff_t,
          class P = V*,
          class R = V&>
struct random_access_iterator_helper
  : random_access_iteratable<T, P, D, R
  , boost::iterator<std::random_access_iterator_tag, V, D, P, R
    > >
{
  friend D requires_difference_operator(const T& x, const T& y) {
    return x - y;
  }
}; // random_access_iterator_helper

*/

} // namespace seqan

#endif // SEQAN_MATH_OPERATORS_H_
