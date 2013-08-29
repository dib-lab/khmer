# Convenience functions for performing common argument-checking tasks in 
# scripts.


def print_error( msg ):
    """
        Prints the given message to 'stderr'.
    """

    import sys

    print >>sys.stderr, msg
    

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
