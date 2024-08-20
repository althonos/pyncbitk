from ...objects.seqalign.seq_align_set cimport CSeq_align_set


cdef extern from "objtools/align_format/align_format_util.hpp" namespace "ncbi::align_format" nogil:

    cppclass CAlignFormatUtil:
        @staticmethod
        int GetMasterCoverage(const CSeq_align_set& alnset)