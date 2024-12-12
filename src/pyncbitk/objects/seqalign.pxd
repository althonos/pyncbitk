# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seqalign.score cimport CScore
from ..toolkit.objects.seqalign.seq_align cimport CSeq_align, TDim, C_Segs
from ..toolkit.objects.seqalign.dense_seg cimport CDense_seg
from ..toolkit.objects.seqalign.seq_align_set cimport CSeq_align_set

from ..serial cimport Serial


cdef class SeqAlignScore(Serial):
    cdef CRef[CScore] _ref

    @staticmethod
    cdef SeqAlignScore _wrap(CRef[CScore] ref)

cdef class AlignRow:
    cdef CRef[CSeq_align] _ref
    cdef TDim             _row

cdef class AlignSegments(Serial):
    cdef CRef[C_Segs]     _ref

    @staticmethod
    cdef AlignSegments _wrap(CRef[C_Segs] ref)

cdef class DenseSegments(AlignSegments):
    pass

cdef class DenseSegmentsData(Serial):
    cdef CRef[CDense_seg] _ref

    @staticmethod
    cdef DenseSegmentsData _wrap(CRef[CDense_seg] ref)


cdef class SeqAlign(Serial):
    cdef CRef[CSeq_align] _ref

    @staticmethod
    cdef SeqAlign _wrap(CRef[CSeq_align] ref)

cdef class GlobalSeqAlign(SeqAlign):
    pass

cdef class DiagonalSeqAlign(SeqAlign):
    pass

cdef class PartialSeqAlign(SeqAlign):
    pass

cdef class DiscontinuousSeqAlign(SeqAlign):
    pass

cdef class SeqAlignSet(Serial):
    cdef CRef[CSeq_align_set] _ref

    @staticmethod
    cdef SeqAlignSet _wrap(CRef[CSeq_align_set] ref)

    cpdef int master_coverage(self) except? 0