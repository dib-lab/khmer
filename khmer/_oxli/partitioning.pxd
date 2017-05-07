from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint8_t, uint64_t
from libc.stdio cimport FILE

from wrapper cimport ComponentPtr, ComponentPtrVector, CpGuardedHashCompMap
from wrapper  cimport CpHashgraph, CpStreamingPartitioner, BoundedCounterType


cdef class Component:
    cdef ComponentPtr _this

    cdef void save(self, FILE * fp)

    @staticmethod
    cdef Component wrap(ComponentPtr ptr)

    @staticmethod
    cdef vector[BoundedCounterType] _tag_counts(ComponentPtr comp, CpHashgraph* graph)

    @staticmethod
    cdef float _mean_tag_count(ComponentPtr comp, CpHashgraph * graph)

    @staticmethod
    cdef ComponentPtr load(uint64_t component_id, list tags)


cdef class StreamingPartitioner:
    cdef shared_ptr[CpStreamingPartitioner] _this
    cdef weak_ptr[ComponentPtrVector] _components
    cdef weak_ptr[CpGuardedHashCompMap] _tag_component_map
    cdef CpHashgraph * _graph_ptr
    cdef readonly object graph
    cdef readonly uint64_t n_consumed

