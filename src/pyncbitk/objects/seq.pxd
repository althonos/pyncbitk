# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.bioseq cimport CBioseq

from ..serial cimport Serial

cdef class BioSeq(Serial):
    cdef CRef[CBioseq] _ref

    @staticmethod
    cdef BioSeq _wrap(CRef[CBioseq] ref)