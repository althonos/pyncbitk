from .ncbitype cimport Int4, Uint4

cdef extern from "corelib/ncbimisc.hpp" namespace "ncbi" nogil:

    ctypedef unsigned int TSeqPos
    const TSeqPos kInvalidSeqPos
    ctypedef int TSignedSeqPos

    ctypedef int TGi
    ctypedef Int4 TIntId
    ctypedef Uint4 TUintId

    cdef struct SStrictId_Gi:
        TIntId TId

    cdef struct SStrictId_Entrez:
        TIntId TId

    cdef struct SStrictId_Tax:
        int TId

    ctypedef int TTaxId