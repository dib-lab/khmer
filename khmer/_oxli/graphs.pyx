# cython: c_string_type=unicode, c_string_encoding=utf8
from cython.operator cimport dereference as deref

from .._khmer import Countgraph, Nodegraph, GraphLabels 


cdef CpHashgraph * get_hashgraph_ptr(object graph):
    if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
        return NULL
    
    cdef CPyHashgraph_Object* ptr = <CPyHashgraph_Object*> graph
    return deref(ptr).hashgraph
    

cdef CpLabelHash * get_labelhash_ptr(object labels):
    if not isinstance(labels, GraphLabels):
        return NULL
    
    cdef CPyGraphLabels_Object * ptr = <CPyGraphLabels_Object*> labels
    return deref(ptr).labelhash
