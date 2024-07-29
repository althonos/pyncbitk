from libcpp cimport bool
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t 

from ....corelib.ncbiobj cimport CObject
from .sseqloc cimport SSeqLoc, TSeqLocVector
from .blast_types cimport EProgram, TSeqAlignVector
from .blast_options cimport CBlastOptions


cdef extern from "algo/blast/api/bl2seq.hpp" namespace "ncbi::blast" nogil:

    cppclass CBl2Seq(CObject):
        CBl2Seq(const SSeqLoc& query, const SSeqLoc& subject, EProgram p)

        CBl2Seq(const SSeqLoc& query, const TSeqLocVector& subjects, EProgram p)
        CBl2Seq(const SSeqLoc& query, const TSeqLocVector& subjects, EProgram p, bool dbscan_mode)
        
        # CBl2Seq(const SSeqLoc& query, const SSeqLoc& subject, CBlastOptionsHandle& opts)
        # CBl2Seq(const SSeqLoc& query, const TSeqLocVector& subject, CBlastOptionsHandle& opts)
        # CBl2Seq(const SSeqLoc& query, const TSeqLocVector& subject, CBlastOptionsHandle& opts, bool dbscan_mode)

        TSeqAlignVector Run() except +