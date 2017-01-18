# cython: c_string_type=unicode, c_string_encoding=utf8
import cython
from cython.operator cimport dereference as deref

from .._khmer import Countgraph
from .._khmer import Nodegraph


cdef CpHashtable * get_hashtable_ptr(object graph):
    if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
        return NULL
    
    cdef CPyHashtable_Object* ptr = <CPyHashtable_Object*> graph
    return deref(ptr).hashtable
    
