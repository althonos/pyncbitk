# cython: language_level=3

from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_inst cimport CSeq_inst
from ..toolkit.objects.seq.delta_seq cimport CDelta_seq

from ..serial cimport Serial


# --- SeqInst ------------------------------------------------------------------

cdef class SeqInst(Serial):
    cdef CRef[CSeq_inst] _ref

    @staticmethod
    cdef SeqInst _wrap(CRef[CSeq_inst] ref)

cdef class VirtualInst(SeqInst):
    pass

cdef class ContinuousInst(SeqInst):
    pass

cdef class SegmentedInst(SeqInst):
    pass

cdef class ConstructedInst(SeqInst):
    pass

cdef class RefInst(SeqInst):
    pass

cdef class ConsensusInst(SeqInst):
    pass

cdef class MapInst(SeqInst):
    pass

cdef class DeltaInst(SeqInst):
    
    cpdef ContinuousInst to_continuous(self)


# --- Delta --------------------------------------------------------------------

cdef class Delta(Serial):
    cdef CRef[CDelta_seq] _ref

    @staticmethod
    cdef Delta _wrap(CRef[CDelta_seq] ref)


cdef class LiteralDelta(Delta):
    pass


cdef class LocDelta(Delta):
    pass