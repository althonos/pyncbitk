# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool

from ..toolkit.objects.seqloc.seq_id cimport CSeq_id
from ..toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from ..toolkit.corelib.ncbiobj cimport CRef

from ..serial cimport Serial

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

cdef class OtherId(SeqId):
    pass

# --- TextSeqId ----------------------------------------------------------------

cdef class TextSeqId(Serial):
    cdef CRef[CTextseq_id] _ref

    @staticmethod
    cdef TextSeqId _wrap(CRef[CTextseq_id] ref)