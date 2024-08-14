# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool

from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id
from ..toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from ..toolkit.corelib.ncbiobj cimport CRef

from ..serial cimport Serial


# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc(Serial):
    cdef CRef[CSeq_loc] _loc

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


# --- SeqId --------------------------------------------------------------------

cdef class SeqId(Serial):
    cdef CRef[CSeq_id] _ref

    @staticmethod
    cdef SeqId _wrap(CRef[CSeq_id] ref)


cdef class LocalId(SeqId):
    pass
        
cdef class RefSeqId(SeqId):
    pass

cdef class GenBankId(SeqId):
    pass
    
cdef class ProteinDataBankId(SeqId):
    pass

cdef class GeneralId(SeqId):
    pass

cdef class TextSeqId(Serial):
    cdef CRef[CTextseq_id] _ref

    @staticmethod
    cdef TextSeqId _wrap(CRef[CTextseq_id] ref)
