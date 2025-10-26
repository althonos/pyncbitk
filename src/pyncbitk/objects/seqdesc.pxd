# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_descr cimport CSeq_descr
from ..toolkit.objects.seq.seqdesc cimport CSeqdesc, E_Choice as CSeqdesc_choice

from ..serial cimport Serial


# --- SeqDesc ------------------------------------------------------------------

cdef class SeqDesc(Serial):
    cdef CRef[CSeqdesc] _ref

    @staticmethod
    cdef SeqDesc _wrap(CRef[CSeqdesc] ref)


cdef class NameDesc(SeqDesc):
    pass


cdef class TitleDesc(SeqDesc):
    pass


cdef class RegionDesc(SeqDesc):
    pass


# --- SeqDescSet ---------------------------------------------------------------

cdef class SeqDescSet(Serial):
    cdef CRef[CSeq_descr] _ref

    @staticmethod
    cdef SeqDescSet _wrap(CRef[CSeq_descr] ref)
