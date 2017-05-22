# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from __future__ import unicode_literals
from cpython.version cimport PY_MAJOR_VERSION


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
        raise TypeError("Requires a string-like sequence")

    if isinstance(s, unicode):
        s = s.encode('utf-8')
    return s


cdef unicode _ustring(s):
    if type(s) is unicode:
        # fast path for most common case(s)
        return <unicode>s
    elif PY_MAJOR_VERSION < 3 and isinstance(s, bytes):
        # only accept byte strings in Python 2.x, not in Py3
        return (<bytes>s).decode('UTF-8')
    elif isinstance(s, unicode):
        # an evil cast to <unicode> might work here in some(!) cases,
        # depending on what the further processing does.  to be safe,
        # we can always create a copy instead
        return unicode(s)
    else:
        raise TypeError(...)
