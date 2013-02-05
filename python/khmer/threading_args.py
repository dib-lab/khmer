"""
    Provides functions to manage threading arguments for the various scripts.
"""

DEFAULT_N_THREADS = 1

def add_threading_args( parser ):
    
    parser.add_argument(
	'--threads', '-T', dest='n_threads',
        default=DEFAULT_N_THREADS,
        help='Number of simultaneous threads to execute'
    )

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
