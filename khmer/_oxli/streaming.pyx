from cython.operator cimport dereference as deref
from libcpp cimport bool

from _oxli cimport CpSequence


cdef class Sequence:

    def __cinit__(self, create=False):
        if create is True:
            self._this.reset(new CpSequence())
    
    @property
    def name(self):
        return deref(self._this).name

    @property
    def annotations(self):
        return deref(self._this).annotations

    @property
    def sequence(self):
        return deref(self._this).sequence

    @property
    def quality(self):
        return deref(self._this).quality

    @staticmethod
    def create(str sequence, str name, str annotations=None, str quality=None):
        cdef Sequence seq = Sequence(create=True)
        deref(seq._this).sequence = sequence
        deref(seq._this).name = name
        if annotations is not None:
            deref(seq._this).annotations = annotations
        if quality is not None:
            deref(seq._this).quality = quality

cdef class SequencePair:
    
    def __cinit__(self, Sequence first, Sequence second):
        self.first = first
        self.second = second


cdef class ReadBundle:

    def __cinit__(self, *raw_records):
        self.reads = [r for r in raw_records if r]

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def total_length(self):
        return sum([len(r.sequence) for r in self.reads])


cdef class FastxParser:

    def __cinit__(self, str filename):
        self._this.reset(new CpFastxParser(filename.encode()))

    cdef Sequence _next(self):
        cdef Sequence seq = Sequence(create=True)
        if deref(self._this).is_complete():
            return None
        else:
            deref(self._this).imprint_next_read(deref(seq._this))
            return seq

    def __iter__(self):
        cdef Sequence seq = self._next()
        while seq is not None:
            yield seq
            seq = self._next()


cpdef tuple broken_paired_reader(FastxParser fastx_iter, 
                                 int min_length=None,
                                 bool force_single=False, 
                                 bool require_paired=False):
    pass


cdef tuple _split_left_right(str name):
    pass


cdef bool check_is_pair(Sequence first, Sequence second):
    pass


cdef bool check_is_left(str name):
    pass


cdef bool check_is_right(str name):
    pass


