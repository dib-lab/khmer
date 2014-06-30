"""
k-mer count and presence
"""

from galaxy.datatypes.binary import Binary

import os
import logging

log = logging.getLogger(__name__)

class Count( Binary ):

    def __init__( self, **kwd ):
        Binary.__init__( self, **kwd )


class Presence( Binary ):

    def __init__( self, **kwd ):
        Binary.__init__( self, **kwd )

Binary.register_unsniffable_binary_ext("ct")
Binary.register_unsniffable_binary_ext("pt")
