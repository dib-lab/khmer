# -*- coding: UTF-8 -*-
from cython.operator cimport dereference as deref
cimport cython

from khmer._oxli.utils cimport _bstring
from khmer._oxli.graphs cimport Hashtable

cdef class Alphabets:
    
    @staticmethod
    def get(name):
        cdef string alphabet = Alphabets._get(name)
        return alphabet

    @staticmethod
    cdef string _get(str name) except *:
        if name == 'DNA_SIMPLE':
            return DNA_SIMPLE
        elif name == 'DNAN_SIMPLE':
            return DNAN_SIMPLE
        elif name == 'RNA_SIMPLE':
            return RNA_SIMPLE
        elif name == 'RNAN_SIMPLE':
            return RNAN_SIMPLE
        elif name == 'IUPAC_NUCL':
            return IUPAC_NUCL
        elif name == 'IUPAC_AA':
            return IUPAC_AA
        else:
            raise ValueError('No alphabet with name {0}'.format(name))


@cython.freelist(100)
cdef class Sequence:

    def __cinit__(self, name=None, sequence=None,
                        quality=None, description=None,
                        cleaned_seq=None):

        if name is not None and sequence is not None:
            self._obj.sequence = _bstring(sequence)
            self._obj.name = _bstring(name)
            if description is not None:
                self._obj.description = _bstring(description)
            if quality is not None:
                self._obj.quality = _bstring(quality)
            if cleaned_seq is not None:
                self._obj.cleaned_seq = _bstring(cleaned_seq)
            else:
                self._obj.cleaned_seq = self._obj.sequence

    def __str__(self):
        return self.cleaned_seq if self._obj.cleaned_seq.length() > 0 else self.sequence

    def __repr__(self):
        return 'Sequence(name="{0}", sequence="{1}")'.format(self.name, self.sequence)

    def __len__(self):
        return self._obj.sequence.length()

    def __richcmp__(x, y, op):
        if op == 2:
            return x.name == y.name and x.sequence == y.sequence
        else:
            raise NotImplementedError('Operator not available')

    def kmers(self, int K):
        cdef int i = 0
        cdef unicode sequence = self.sequence
        for i in range(0, len(self)-K+1):
            yield sequence[i:i+K]

    def __getitem__(self, x):
        # Definitely optimize this.
        return self.sequence[x]

    def trim(self, int trim_at):
        self._obj.sequence.resize(trim_at)
        self._obj.cleaned_seq.resize(trim_at)
        if self._obj.quality.length() != 0:
            self._obj.quality.resize(trim_at)

    def clean(self):
        '''Calls set_cleaned_seq() on the underlying container.'''
        self._obj.set_clean_seq()

    @property
    def name(self):
        cdef unicode name = self._obj.name
        return name if name else None

    @property
    def sequence(self):
        cdef unicode sequence = self._obj.sequence
        return sequence if sequence else None

    @property
    def description(self):
        cdef unicode description = self._obj.description
        return description if description else None

    @property
    def quality(self):
        cdef unicode quality = self._obj.quality
        return quality if quality else None

    @property
    def cleaned_seq(self):
        cdef unicode cleaned_seq = self._obj.cleaned_seq
        return cleaned_seq if cleaned_seq else None

    @staticmethod
    def from_screed_record(record):
        cdef Sequence seq = Sequence(name=record.name,
                                     sequence=record.sequence)
        if hasattr(record, 'quality'):
            seq._obj.quality = _bstring(record.quality)

        for attr in ('annotations', 'description'):
            if hasattr(record, attr):
                seq._obj.description = _bstring(getattr(record, attr))

        return seq

    @staticmethod
    cdef Sequence _wrap(CpSequence cseq):
        cdef Sequence seq = Sequence()
        seq._obj = cseq
        return seq


cdef string _object_to_string(object sequence) except *:
    if isinstance(sequence, bytes):
        return sequence
    elif isinstance(sequence, Sequence):
        return (<Sequence>sequence)._obj.cleaned_seq
    else:
        return _bstring(sequence)


cdef class ReadBundle:

    def __cinit__(self, *raw_records):
        self.reads = [r for r in raw_records if r]

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def total_length(self):
        return sum([len(r.sequence) for r in self.reads])


cdef bool is_valid(const char base, string& alphabet):
    cdef char b
    for b in alphabet:
        if b == base:
            return True
    return False


cdef bool sanitize_sequence(string& sequence,
                                   string& alphabet,
                                   bool convert_n):
    cdef int i = 0
    for i in range(sequence.length()):
        sequence[i] &= 0xdf
        if not is_valid(sequence[i], alphabet):
            return False
        if convert_n and sequence[i] == b'N':
            sequence[i] = b'A'
    return True


def trim_sequence(Hashtable graph, Sequence record, int cutoff,
                  variable_coverage=False, normalize_to=None):
    if variable_coverage:
        if not graph.median_at_least(record.cleaned_seq, normalize_to):
            return record, False

    trim_at = graph._trim_on_abundance(record, cutoff)
    
    if trim_at < graph.ksize():
        return None, True

    if trim_at == len(record):
        return record, False

    record.trim(trim_at)
    return record, True
