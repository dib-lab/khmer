
def is_prime(n):
    return _is_prime(n)


def get_n_primes_near_x(n_primes, x):
    primes = _get_n_primes_near_x(n_primes, x)
    if len(primes) != n_primes:
        msg = "unable to find {0} prime numbers < {1}".format(n_primes, x)
        raise RuntimeError(msg)
    return primes 
