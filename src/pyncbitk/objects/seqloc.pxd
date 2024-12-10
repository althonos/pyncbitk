# cython: language_level=3

from libcpp cimport bool

from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc
from ..toolkit.corelib.ncbiobj cimport CRef

from ..serial cimport Serial


# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc(Serial):
    cdef CRef[CSeq_loc] _loc

    @staticmethod
    cdef SeqLoc _wrap(CRef[CSeq_loc] ref)

cdef class NullLoc(SeqLoc):
    pass

cdef class EmptySeqLoc(SeqLoc):
    pass

cdef class WholeSeqLoc(SeqLoc):
    pass

cdef class SeqIntervalLoc(SeqLoc):
    pass

cdef class PackedSeqLoc(SeqLoc):
    pass

cdef class PackedPointsLoc(SeqLoc):
    pass

cdef class MixLoc(SeqLoc):
    pass

cdef class EquivalentLoc(SeqLoc):
    pass

cdef class BondLoc(SeqLoc):
    pass

cdef class FeatureLoc(SeqLoc):
    pass
