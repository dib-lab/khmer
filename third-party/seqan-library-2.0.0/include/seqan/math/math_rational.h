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
// A data type to represent Rational numbers.
// Taken from the Rational Number Library in Boost version 1.47.
// ==========================================================================

#ifndef SEQAN_MATH_RATIONAL_H_
#define SEQAN_MATH_RATIONAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TInt>
class Rational;

template <typename TInt>
Rational<TInt> abs(const Rational<TInt>& r);


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TInt>
class Rational :
    less_than_comparable < Rational<TInt>,
    equality_comparable < Rational<TInt>,
    less_than_comparable2 < Rational<TInt>, TInt,
    equality_comparable2 < Rational<TInt>, TInt,
    addable < Rational<TInt>,
    subtractable < Rational<TInt>,
    multipliable < Rational<TInt>,
    dividable < Rational<TInt>,
    addable2 < Rational<TInt>, TInt,
    subtractable2 < Rational<TInt>, TInt,
    subtractable2_left < Rational<TInt>, TInt,
    multipliable2 < Rational<TInt>, TInt,
    dividable2 < Rational<TInt>, TInt,
    dividable2_left < Rational<TInt>, TInt,
    incrementable < Rational<TInt>,
    decrementable < Rational<TInt>
    > > > > > > > > > > > > > > > >
{
    // Class-wide pre-conditions
//    SEQAN_STATIC_ASSERT( std::numeric_limits<TInt>::is_specialized );

    // Helper types
//    typedef typename boost::call_traits<TInt>::param_type param_type;
    typedef const TInt & param_type;

    struct helper { TInt parts[2]; };
    typedef TInt (helper::* bool_type)[2];

public:
    typedef TInt int_type;
    Rational() : num(0), den(1) {}

    template <typename T>
    Rational(T const & n, SEQAN_CTOR_ENABLE_IF( IsInteger<T> ) ) : num(n), den(1) { (void)dummy; }
    Rational(param_type n, param_type d) : num(n), den(d) { normalize(); }

    // Default copy constructor and assignment are fine

    // Add assignment from TInt
    Rational& operator=(param_type n) { return assign(n, 1); }

    // Assign in place
    Rational& assign(param_type n, param_type d);

    // Access to representation
    TInt numerator() const { return num; }
    TInt denominator() const { return den; }

    // Arithmetic assignment operators
    Rational& operator+= (const Rational& r);
    Rational& operator-= (const Rational& r);
    Rational& operator*= (const Rational& r);
    Rational& operator/= (const Rational& r);

    Rational& operator+= (param_type i);
    Rational& operator-= (param_type i);
    Rational& operator*= (param_type i);
    Rational& operator/= (param_type i);

    // Increment and decrement
    const Rational& operator++();
    const Rational& operator--();

    // Operator not
    bool operator!() const { return !num; }

    // Boolean conversion


    operator bool_type() const { return operator !() ? 0 : &helper::parts; }
    operator double() const
    {
        SEQAN_ASSERT_NEQ (den, TInt(0));
        return (double)num / (double)den;
    }

    // Comparison operators
    bool operator< (const Rational& r) const;
    bool operator== (const Rational& r) const;

    bool operator< (param_type i) const;
    bool operator> (param_type i) const;
    bool operator== (param_type i) const;

private:
    // Implementation - numerator and denominator (normalized).
    // Other possibilities - separate whole-part, or sign, fields?
    TInt num;
    TInt den;

    // Representation note: Fractions are kept in normalized form at all
    // times. normalized form is defined as gcd(num,den) == 1 and den > 0.
    // In particular, note that the implementation of abs() below relies
    // on den always being positive.
    bool test_invariant() const;
    void normalize();
};

// Assign in place
template <typename TInt>
inline Rational<TInt>& Rational<TInt>::assign(param_type n, param_type d)
{
    num = n;
    den = d;
    normalize();
    return *this;
}

// Unary plus and minus
template <typename TInt>
inline Rational<TInt> operator+ (const Rational<TInt>& r)
{
    return r;
}

template <typename TInt>
inline Rational<TInt> operator- (const Rational<TInt>& r)
{
    return Rational<TInt>(-r.numerator(), r.denominator());
}

// Arithmetic assignment operators
template <typename TInt>
Rational<TInt>& Rational<TInt>::operator+= (const Rational<TInt>& r)
{
    // This calculation avoids overflow, and minimises the number of expensive
    // calculations. Thanks to Nickolay Mladenov for this algorithm.
    //
    // Proof:
    // We have to compute a/b + c/d, where gcd(a,b)=1 and gcd(b,c)=1.
    // Let g = gcd(b,d), and b = b1*g, d=d1*g. Then gcd(b1,d1)=1
    //
    // The result is (a*d1 + c*b1) / (b1*d1*g).
    // Now we have to normalize this ratio.
    // Let's assume h | gcd((a*d1 + c*b1), (b1*d1*g)), and h > 1
    // If h | b1 then gcd(h,d1)=1 and hence h|(a*d1+c*b1) => h|a.
    // But since gcd(a,b1)=1 we have h=1.
    // Similarly h|d1 leads to h=1.
    // So we have that h | gcd((a*d1 + c*b1) , (b1*d1*g)) => h|g
    // Finally we have gcd((a*d1 + c*b1), (b1*d1*g)) = gcd((a*d1 + c*b1), g)
    // Which proves that instead of normalizing the result, it is better to
    // divide num and den by gcd((a*d1 + c*b1), g)

    // Protect against self-modification
    TInt r_num = r.num;
    TInt r_den = r.den;

    TInt g = greatestCommonDivisor(den, r_den);
    den /= g;  // = b1 from the calculations above
    num = num * (r_den / g) + r_num * den;
    g = greatestCommonDivisor(num, g);
    num /= g;
    den *= r_den/g;

    return *this;
}

template <typename TInt>
Rational<TInt>& Rational<TInt>::operator-= (const Rational<TInt>& r)
{
    // Protect against self-modification
    TInt r_num = r.num;
    TInt r_den = r.den;

    // This calculation avoids overflow, and minimises the number of expensive
    // calculations. It corresponds exactly to the += case above
    TInt g = greatestCommonDivisor(den, r_den);
    den /= g;
    num = num * (r_den / g) - r_num * den;
    g = greatestCommonDivisor(num, g);
    num /= g;
    den *= r_den/g;

    return *this;
}

template <typename TInt>
Rational<TInt>& Rational<TInt>::operator*= (const Rational<TInt>& r)
{
    // Protect against self-modification
    TInt r_num = r.num;
    TInt r_den = r.den;

    // Avoid overflow and preserve normalization
    TInt gcd1 = greatestCommonDivisor(num, r_den);
    TInt gcd2 = greatestCommonDivisor(r_num, den);
    num = (num/gcd1) * (r_num/gcd2);
    den = (den/gcd2) * (r_den/gcd1);
    return *this;
}

template <typename TInt>
Rational<TInt>& Rational<TInt>::operator/= (const Rational<TInt>& r)
{
    // Protect against self-modification
    TInt r_num = r.num;
    TInt r_den = r.den;

    // Avoid repeated construction
    TInt zero(0);

    // Trap division by zero
    SEQAN_ASSERT_NEQ (r_num, zero);

    if (num == zero)
        return *this;

    // Avoid overflow and preserve normalization
    TInt gcd1 = greatestCommonDivisor(num, r_num);
    TInt gcd2 = greatestCommonDivisor(r_den, den);
    num = (num/gcd1) * (r_den/gcd2);
    den = (den/gcd2) * (r_num/gcd1);

    if (den < zero) {
        num = -num;
        den = -den;
    }
    return *this;
}

// Mixed-mode operators
template <typename TInt>
inline Rational<TInt>&
Rational<TInt>::operator+= (param_type i)
{
    return operator+= (Rational<TInt>(i));
}

template <typename TInt>
inline Rational<TInt>&
Rational<TInt>::operator-= (param_type i)
{
    return operator-= (Rational<TInt>(i));
}

template <typename TInt>
inline Rational<TInt>&
Rational<TInt>::operator*= (param_type i)
{
    return operator*= (Rational<TInt>(i));
}

template <typename TInt>
inline Rational<TInt>&
Rational<TInt>::operator/= (param_type i)
{
    return operator/= (Rational<TInt>(i));
}

// Increment and decrement
template <typename TInt>
inline const Rational<TInt>& Rational<TInt>::operator++()
{
    // This can never denormalise the fraction
    num += den;
    return *this;
}

template <typename TInt>
inline const Rational<TInt>& Rational<TInt>::operator--()
{
    // This can never denormalise the fraction
    num -= den;
    return *this;
}

// Comparison operators
template <typename TInt>
bool Rational<TInt>::operator< (const Rational<TInt>& r) const
{
    // Avoid repeated construction
    int_type const  zero( 0 );

    // This should really be a class-wide invariant.  The reason for these
    // checks is that for 2's complement systems, INT_MIN has no corresponding
    // positive, so negating it during normalization keeps it INT_MIN, which
    // is bad for later calculations that assume a positive denominator.
    SEQAN_ASSERT_GT( this->den, zero );
    SEQAN_ASSERT_GT( r.den, zero );

    // Determine relative order by expanding each value to its simple continued
    // fraction representation using the Euclidian GCD algorithm.
    struct { int_type  n, d, q, r; }  ts = { this->num, this->den, this->num /
     this->den, this->num % this->den }, rs = { r.num, r.den, r.num / r.den,
     r.num % r.den };
    unsigned  reverse = 0u;

    // Normalize negative moduli by repeatedly adding the (positive) denominator
    // and decrementing the quotient.  Later cycles should have all positive
    // values, so this only has to be done for the first cycle.  (The rules of
    // C++ require a nonnegative quotient & remainder for a nonnegative dividend
    // & positive divisor.)
    while ( ts.r < zero )  { ts.r += ts.d; --ts.q; }
    while ( rs.r < zero )  { rs.r += rs.d; --rs.q; }

    // Loop through and compare each variable's continued-fraction components
    while ( true )
    {
        // The quotients of the current cycle are the continued-fraction
        // components.  Comparing two c.f. is comparing their sequences,
        // stopping at the first difference.
        if ( ts.q != rs.q )
        {
            // Since reciprocation changes the relative order of two variables,
            // and c.f. use reciprocals, the less/greater-than test reverses
            // after each index.  (Start w/ non-reversed @ whole-number place.)
            return reverse ? ts.q > rs.q : ts.q < rs.q;
        }

        // Prepare the next cycle
        reverse ^= 1u;

        if ( (ts.r == zero) || (rs.r == zero) )
        {
            // At least one variable's c.f. expansion has ended
            break;
        }

        ts.n = ts.d;         ts.d = ts.r;
        ts.q = ts.n / ts.d;  ts.r = ts.n % ts.d;
        rs.n = rs.d;         rs.d = rs.r;
        rs.q = rs.n / rs.d;  rs.r = rs.n % rs.d;
    }

    // Compare infinity-valued components for otherwise equal sequences
    if ( ts.r == rs.r )
    {
        // Both remainders are zero, so the next (and subsequent) c.f.
        // components for both sequences are infinity.  Therefore, the sequences
        // and their corresponding values are equal.
        return false;
    }
    else
    {
        // Exactly one of the remainders is zero, so all following c.f.
        // components of that variable are infinity, while the other variable
        // has a finite next c.f. component.  So that other variable has the
        // lesser value (modulo the reversal flag!).
        return ( ts.r != zero ) != static_cast<bool>( reverse );
    }
}

template <typename TInt>
bool Rational<TInt>::operator< (param_type i) const
{
    // Avoid repeated construction
    int_type const  zero( 0 );

    // Break value into mixed-fraction form, w/ always-nonnegative remainder
    SEQAN_ASSERT_GT( this->den, zero );
    int_type  q = this->num / this->den, r = this->num % this->den;
    while ( r < zero )  { r += this->den; --q; }

    // Compare with just the quotient, since the remainder always bumps the
    // value up.  [Since q = floor(n/d), and if n/d < i then q < i, if n/d == i
    // then q == i, if n/d == i + r/d then q == i, and if n/d >= i + 1 then
    // q >= i + 1 > i; therefore n/d < i iff q < i.]
    return q < i;
}

template <typename TInt>
bool Rational<TInt>::operator> (param_type i) const
{
    // Trap equality first
    if (num == i && den == TInt(1))
        return false;

    // Otherwise, we can use operator<
    return !operator<(i);
}

template <typename TInt>
inline bool Rational<TInt>::operator== (const Rational<TInt>& r) const
{
    return ((num == r.num) && (den == r.den));
}

template <typename TInt>
inline bool Rational<TInt>::operator== (param_type i) const
{
    return ((den == TInt(1)) && (num == i));
}

// Invariant check
template <typename TInt>
inline bool Rational<TInt>::test_invariant() const
{
    return ( this->den > int_type(0) ) && ( greatestCommonDivisor(this->num, this->den) ==
     int_type(1) );
}

// Normalisation
template <typename TInt>
void Rational<TInt>::normalize()
{
    // Avoid repeated construction
    TInt zero(0);

    SEQAN_ASSERT_NEQ (den, zero);

    // Handle the case of zero separately, to avoid division by zero
    if (num == zero) {
        den = TInt(1);
        return;
    }

    TInt g = greatestCommonDivisor(num, den);

    num /= g;
    den /= g;

    // Ensure that the denominator is positive
    if (den < zero) {
        num = -num;
        den = -den;
    }

    SEQAN_ASSERT( this->test_invariant() );
}

// Input and output
template <typename TInt>
std::istream& operator>> (std::istream& is, Rational<TInt>& r)
{
    TInt n = TInt(0), d = TInt(1);
    char c = 0;

    is >> n;
    if (!is) return is;

    c = is.get();
    if (c == '/')
    {
        is >> std::noskipws >> d;
    }
    else if (c == '.')
    {
        bool negative = false;
        if (n < TInt(0))
        {
            n = -n;
            negative = true;
        }
        c = is.get();
        // read digits as long we can store them
        while ('0' <= c && c <= '9' &&
                (n < (TInt)MaxValue<TInt>::VALUE / (TInt)10 - (TInt)9) &&
                (d < (TInt)MaxValue<TInt>::VALUE / (TInt)10))
        {
            n = 10 * n + (c - '0');
            d *= 10;
            c = is.get();
        }
        // ignore remaining digits
        while ('0' <= c && c <= '9')
            is.get();
        is.unget();

        if (negative) n = -n;
    }
    r.assign(n, d);

    return is;
}

// Add manipulators for output format?
template <typename TInt>
std::ostream& operator<< (std::ostream& os, const Rational<TInt>& r)
{
    os << r.numerator() << '/' << r.denominator();
    return os;
}

// Type conversion
template <typename T, typename TInt>
inline T rational_cast(const Rational<TInt>& src)
{
    return static_cast<T>(src.numerator())/static_cast<T>(src.denominator());
}

// Do not use any abs() defined on TInt - it isn't worth it, given the
// difficulties involved (Koenig lookup required, there may not *be* an abs()
// defined, etc etc).
template <typename TInt>
inline Rational<TInt> abs(const Rational<TInt>& r)
{
    if (r.numerator() >= TInt(0))
        return r;

    return Rational<TInt>(-r.numerator(), r.denominator());
}


template <typename TInt>
inline TInt floor(const Rational<TInt>& r)
{
    // Avoid repeated construction
    TInt zero(0);

    SEQAN_ASSERT_NEQ (r.denominator(), zero);

    if (r.numerator() >= zero)
        return r.numerator() / r.denominator();
    else
        return ((r.numerator() + 1) / r.denominator()) - 1;
}

template <typename TInt>
inline TInt ceil(const Rational<TInt>& r)
{
    // Avoid repeated construction
    TInt zero(0);

    SEQAN_ASSERT_NEQ (r.denominator(), zero);

    if (r.numerator() > zero)
        return ((r.numerator() - 1) / r.denominator()) + 1;
    else
        return r.numerator() / r.denominator();
}


}  // namespace seqan

#endif  // #ifndef SEQAN_MATH_RATIONAL_H_
