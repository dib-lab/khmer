import cython
from cython.operator cimport dereference as deref, preincrement as inc
from libc.limits cimport UINT_MAX

cdef class Component:

    cdef ComponentPtr _this

    def __cinit__(self, Component other=None):
        if other is not None:
            self._this.reset(other._this.get())
        else:
            self._this.reset(new CyComponent())
        
    property component_id:
        def __get__(self):
            return deref(self._this).component_id

    property n_merges:
        def __get__(self):
            return deref(self._this).get_n_merges()

    def __len__(self):
        return deref(self._this).get_n_tags()

    def __iter__(self):
        it = deref(self._this).tags.begin()
        while it != deref(self._this).tags.end():
            yield deref(it)
            inc(it)

