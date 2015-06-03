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
// A data type to represent rational numbers.
// Taken from the Rational Number Library in Boost version 1.47 adapted
// to SeqAn's code conventions.
// ==========================================================================

#ifndef SEQAN_MATH_COMMON_FACTOR_H_
#define SEQAN_MATH_COMMON_FACTOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template < typename IntegerType >
inline IntegerType
greatestCommonDivisor( IntegerType const &a, IntegerType const &b );

template < typename IntegerType >
inline IntegerType
leastCommonMultiple( IntegerType const &a, IntegerType const &b );


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template < typename IntegerType >
class GcdEvaluator
{
public:
    // Types
    typedef IntegerType  resultType, firstArgumentType, secondArgumentType;

    // Function object interface
    resultType  operator ()( firstArgumentType const &a,
     secondArgumentType const &b ) const;

};  // boost::math::GcdEvaluator


//  Least common multiple evaluator class declaration  -----------------------//

template < typename IntegerType >
class LcmEvaluator
{
public:
    // Types
    typedef IntegerType  resultType, firstArgumentType, secondArgumentType;

    // Function object interface
    resultType  operator ()( firstArgumentType const &a,
     secondArgumentType const &b ) const;

};  // boost::math::LcmEvaluator

//  Implementation details  --------------------------------------------------//

namespace detail
{
    // Greatest common divisor for rings (including unsigned integers)
    template < typename RingType >
    RingType
    gcdEuclidean
    (
        RingType a,
        RingType b
    )
    {
        // Avoid repeated construction
        RingType const  zero = static_cast<RingType>( 0 );

        // Reduce by GCD-remainder property [GCD(a,b) == GCD(b,a MOD b)]
        while ( true )
        {
            if ( a == zero )
                return b;
            b %= a;

            if ( b == zero )
                return a;
            a %= b;
        }
    }

    // Greatest common divisor for (signed) integers
    template < typename IntegerType >
    inline
    IntegerType
    gcdInteger
    (
        IntegerType const &  a,
        IntegerType const &  b
    )
    {
        // Avoid repeated construction
        IntegerType const  zero = static_cast<IntegerType>( 0 );
        IntegerType const  result = gcdEuclidean( a, b );

        return ( result < zero ) ? static_cast<IntegerType>(-result) : result;
    }

    // Greatest common divisor for unsigned binary integers
    template < typename BuiltInUnsigned >
    BuiltInUnsigned
    gcd_binary
    (
        BuiltInUnsigned  u,
        BuiltInUnsigned  v
    )
    {
        if ( u && v )
        {
            // Shift out common factors of 2
            unsigned  shifts = 0;

            while ( !(u & 1u) && !(v & 1u) )
            {
                ++shifts;
                u >>= 1;
                v >>= 1;
            }

            // Start with the still-even one, if any
            BuiltInUnsigned  r[] = { u, v };
            unsigned         which = static_cast<bool>( u & 1u );

            // Whittle down the values via their differences
            do
            {
                // Remove factors of two from the even one
                while ( !(r[ which ] & 1u) )
                {
                    r[ which ] >>= 1;
                }

                // Replace the larger of the two with their difference
                if ( r[!which] > r[which] )
                {
                    which ^= 1u;
                }

                r[ which ] -= r[ !which ];
            }
            while ( r[which] );

            // Shift-in the common factor of 2 to the residues' GCD
            return r[ !which ] << shifts;
        }
        else
        {
            // At least one input is zero, return the other
            // (adding since zero is the additive identity)
            // or zero if both are zero.
            return u + v;
        }
    }

    // Least common multiple for rings (including unsigned integers)
    template < typename RingType >
    inline
    RingType
    lcmEuclidean
    (
        RingType const &  a,
        RingType const &  b
    )
    {
        RingType const  zero = static_cast<RingType>( 0 );
        RingType const  temp = gcdEuclidean( a, b );

        return ( temp != zero ) ? ( a / temp * b ) : zero;
    }

    // Least common multiple for (signed) integers
    template < typename IntegerType >
    inline
    IntegerType
    lcmInteger
    (
        IntegerType const &  a,
        IntegerType const &  b
    )
    {
        // Avoid repeated construction
        IntegerType const  zero = static_cast<IntegerType>( 0 );
        IntegerType const  result = lcmEuclidean( a, b );

        return ( result < zero ) ? static_cast<IntegerType>(-result) : result;
    }

    // Function objects to find the best way of computing GCD or LCM
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    template < typename T, bool IsSpecialized, bool IsSigned >
    struct GcdOptimalEvaluatorHelper
    {
        T  operator ()( T const &a, T const &b )
        {
            return gcdEuclidean( a, b );
        }
    };

    template < typename T >
    struct GcdOptimalEvaluatorHelper< T, true, true >
    {
        T  operator ()( T const &a, T const &b )
        {
            return gcdInteger( a, b );
        }
    };
#else
    template < bool IsSpecialized, bool IsSigned >
    struct GcdOptimalEvaluatorHelper2
    {
        template < typename T >
        struct Helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return gcdEuclidean( a, b );
            }
        };
    };

    template < >
    struct GcdOptimalEvaluatorHelper2< true, true >
    {
        template < typename T >
        struct Helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return gcdInteger( a, b );
            }
        };
    };

    template < typename T, bool IsSpecialized, bool IsSigned >
    struct GcdOptimalEvaluatorHelper
        : GcdOptimalEvaluatorHelper2<IsSpecialized, IsSigned>
           ::BOOST_NESTED_TEMPLATE helper<T>
    {
    };
#endif

    template < typename T >
    struct GcdOptimalEvaluator
    {
        T  operator ()( T const &a, T const &b )
        {
            typedef std::numeric_limits<T>  limitsType;

            typedef GcdOptimalEvaluatorHelper<T,
             limitsType::is_specialized, limitsType::is_signed>  helperType;

            helperType  solver;

            return solver( a, b );
        }
    };
#else // BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    template < typename T >
    struct GcdOptimalEvaluator
    {
        T  operator ()( T const &a, T const &b )
        {
            return gcdInteger( a, b );
        }
    };
#endif

    // Specialize for the built-in integers
#define BOOST_PRIVATE_GCD_UF( Ut )                  \
    template < >  struct GcdOptimalEvaluator<Ut>  \
    {  Ut  operator ()( Ut a, Ut b ) const  { return gcd_binary( a, b ); }  }

    BOOST_PRIVATE_GCD_UF( unsigned char );
    BOOST_PRIVATE_GCD_UF( unsigned short );
    BOOST_PRIVATE_GCD_UF( unsigned );
    BOOST_PRIVATE_GCD_UF( unsigned long );

#ifdef BOOST_HAS_LONG_LONG
    BOOST_PRIVATE_GCD_UF( boost::ulong_longType );
#elif defined(BOOST_HAS_MS_INT64)
    BOOST_PRIVATE_GCD_UF( unsigned __int64 );
#endif

#if CHAR_MIN == 0
    BOOST_PRIVATE_GCD_UF( char ); // char is unsigned
#endif

#undef BOOST_PRIVATE_GCD_UF

#define BOOST_PRIVATE_GCD_SF( St, Ut )                            \
    template < >  struct GcdOptimalEvaluator<St>                \
    {  St  operator ()( St a, St b ) const  { Ut const  a_abs =   \
    static_cast<Ut>( a < 0 ? -a : +a ), b_abs = static_cast<Ut>(  \
    b < 0 ? -b : +b ); return static_cast<St>(                    \
    GcdOptimalEvaluator<Ut>()(a_abs, b_abs) ); }  }

    BOOST_PRIVATE_GCD_SF( signed char, unsigned char );
    BOOST_PRIVATE_GCD_SF( short, unsigned short );
    BOOST_PRIVATE_GCD_SF( int, unsigned );
    BOOST_PRIVATE_GCD_SF( long, unsigned long );

#if CHAR_MIN < 0
    BOOST_PRIVATE_GCD_SF( char, unsigned char ); // char is signed
#endif

#ifdef BOOST_HAS_LONG_LONG
    BOOST_PRIVATE_GCD_SF( boost::long_longType, boost::ulong_longType );
#elif defined(BOOST_HAS_MS_INT64)
    BOOST_PRIVATE_GCD_SF( __int64, unsigned __int64 );
#endif

#undef BOOST_PRIVATE_GCD_SF

#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    template < typename T, bool IsSpecialized, bool IsSigned >
    struct LcmOptimalEvaluatorHelper
    {
        T  operator ()( T const &a, T const &b )
        {
            return lcmEuclidean( a, b );
        }
    };

    template < typename T >
    struct LcmOptimalEvaluatorHelper< T, true, true >
    {
        T  operator ()( T const &a, T const &b )
        {
            return lcmInteger( a, b );
        }
    };
#else
    template < bool IsSpecialized, bool IsSigned >
    struct LcmOptimalEvaluatorHelper2
    {
        template < typename T >
        struct Helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return lcmEuclidean( a, b );
            }
        };
    };

    template < >
    struct LcmOptimalEvaluatorHelper2< true, true >
    {
        template < typename T >
        struct Helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return lcmInteger( a, b );
            }
        };
    };

    template < typename T, bool IsSpecialized, bool IsSigned >
    struct LcmOptimalEvaluatorHelper
        : lcmOptimalEvaluatorHelper2<IsSpecialized, IsSigned>
           ::BOOST_NESTED_TEMPLATE helper<T>
    {
    };
#endif

    template < typename T >
    struct LcmOptimalEvaluator
    {
        T  operator ()( T const &a, T const &b )
        {
            typedef std::numeric_limits<T>  limitsType;

            typedef LcmOptimalEvaluatorHelper<T,
             limitsType::is_specialized, limitsType::is_signed>  helperType;

            helperType  solver;

            return solver( a, b );
        }
    };
#else // BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    template < typename T >
    struct LcmOptimalEvaluator
    {
        T  operator ()( T const &a, T const &b )
        {
            return lcmInteger( a, b );
        }
    };
#endif

    // Functions to find the GCD or LCM in the best way
    template < typename T >
    inline
    T
    gcdOptimal
    (
        T const &  a,
        T const &  b
    )
    {
        GcdOptimalEvaluator<T>  solver;

        return solver( a, b );
    }

    template < typename T >
    inline
    T
    lcmOptimal
    (
        T const &  a,
        T const &  b
    )
    {
        LcmOptimalEvaluator<T>  solver;

        return solver( a, b );
    }

}  // namespace detail


//  Greatest common divisor evaluator member function definition  ------------//

template < typename IntegerType >
inline
typename GcdEvaluator<IntegerType>::resultType
GcdEvaluator<IntegerType>::operator ()
(
    firstArgumentType const &   a,
    secondArgumentType const &  b
) const
{
    return detail::gcdOptimal( a, b );
}


//  Least common multiple evaluator member function definition  --------------//

template < typename IntegerType >
inline
typename LcmEvaluator<IntegerType>::resultType
LcmEvaluator<IntegerType>::operator ()
(
    firstArgumentType const &   a,
    secondArgumentType const &  b
) const
{
    return detail::lcmOptimal( a, b );
}


//  Greatest common divisor and least common multiple function definitions  --//

template < typename IntegerType >
inline IntegerType
greatestCommonDivisor
(
    IntegerType const &  a,
    IntegerType const &  b
)
{
    GcdEvaluator<IntegerType>  solver;

    return solver( a, b );
}

template < typename IntegerType >
inline IntegerType
leastCommonMultiple
(
    IntegerType const &  a,
    IntegerType const &  b
)
{
    LcmEvaluator<IntegerType>  solver;

    return solver( a, b );
}

}  // namespace seqan

#endif  // #ifndef SEQAN_MATH_COMMON_FACTOR_H_
