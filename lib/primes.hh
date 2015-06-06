//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef PRIMES_HH
#define PRIMES_HH

#include <math.h>
#include "hashtable.hh"

/* Primes is a class for generating prime numbers. */

namespace khmer
{
class Primes
{

private:
    HashIntoType n;

    /* Returns true if n is prime, false otherwise */
    bool is_prime()
    {
        if (n < 2) {
            return false;
        }
        if (n == 2) {
            return true;
        }
        if ((n % 2) == 0) {
            return false;
        }

        HashIntoType max = (HashIntoType)pow(n, 0.5);
        for (HashIntoType x = 3; x <= max; x+=2)
            if ((n % x) == 0) {
                return false;
            }
        return true;
    }

public:
    Primes(HashIntoType num)
    {
        /* Make sure that the initial number to start from is odd
         * and strictly greater than num */
        n = ++num;
        if ((n % 2) == 0) {
            n++;
        }
    }

    /* Returns the next prime >= n */
    HashIntoType get_next_prime()
    {
        for (;;) {
            if (is_prime()) {
                /* If n is prime, we need to make sure to increment
                 * n by 2 before returning, otherwise we will always
                 * return the same prime number */
                HashIntoType prime = n;
                n += 2;
                return prime;
            }
            n += 2;
        }
    }
};
};

#endif

