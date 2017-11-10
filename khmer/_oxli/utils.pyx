# -*- coding: UTF-8 -*-

from cpython.version cimport PY_MAJOR_VERSION

from cython import short, int, long


def is_prime(n):
    return _is_prime(n)


def get_n_primes_near_x(n_primes, x):
    primes = _get_n_primes_near_x(n_primes, x)
    if len(primes) != n_primes:
        msg = "unable to find {0} prime numbers < {1}".format(n_primes, x)
        raise RuntimeError(msg)
    return primes


cdef bytes _bstring(s):
    if not isinstance(s, (basestring, bytes)):
        raise TypeError("Requires a string-like sequence, "\
                        " got {0} of type {1}".format(s, type(s)))

    if isinstance(s, unicode):
        s = s.encode('utf-8')
    return s


cdef unicode _ustring(s):
    if type(s) is unicode:
        # fast path for most common case(s)
        return <unicode>s
    elif isinstance(s, unicode):
        # an evil cast to <unicode> might work here in some(!) cases,
        # depending on what the further processing does.  to be safe,
        # we can always create a copy instead
        return unicode(s)
    else:
        raise TypeError(...)


cpdef bool is_str(object s):
    return isinstance(s, (basestring, bytes))

cpdef bool is_num(object n):
    return isinstance(n, (int, long))

cdef void _flatten_fill(double * fill_to, object fill_from):
    '''UNSAFE fill from multilevel python iterable to C array.'''
    cdef list flattened = [x for sublist in fill_from for x in sublist]
    for idx, item in enumerate(flattened):
        fill_to[idx] = <double>item

cdef void _fill(double * fill_to, object fill_from):
    '''UNSAFE fill from flat python iterable to C array.'''
    for idx, item in enumerate(fill_from):
        fill_to[idx] = <double>item
