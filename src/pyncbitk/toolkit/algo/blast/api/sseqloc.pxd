from libcpp cimport bool
from libcpp.vector cimport vector

from ....corelib.ncbitype cimport Uint4 
from ....corelib.ncbimisc cimport TSeqPos
from ....corelib.ncbiobj cimport CRef, CConstRef, CObject
from ....objects.seqloc.seq_loc cimport CSeq_loc
from ....objects.seqloc.seq_id cimport CSeq_id
from ....objmgr.scope cimport CScope

cdef extern from "algo/blast/api/sseqloc.hpp" namespace "ncbi::blast" nogil:

    cppclass SSeqLoc:
        CConstRef[CSeq_loc] seqloc
        CRef[CScope] scope
        CRef[CSeq_loc] mask

        bool        ignore_strand_in_mask
        Uint4       genetic_code_id

        SSeqLoc()
        SSeqLoc(const CSeq_loc* sl, CScope* s)
        SSeqLoc(const CSeq_loc* sl, CScope* s, CSeq_loc* m)
        SSeqLoc(const CSeq_loc* sl, CScope* s, CSeq_loc* m, bool ignore_mask_strand)
        SSeqLoc(const CSeq_loc& sl, CScope& s)
        SSeqLoc(const CSeq_loc& sl, CScope& s, CSeq_loc& m)
        SSeqLoc(const CSeq_loc& sl, CScope& s, CSeq_loc& m, bool ignore_mask_strand)

    ctypedef vector[SSeqLoc] TSeqLocVector


    cppclass CBlastSearchQuery(CObject):
        CBlastSearchQuery()
        CBlastSearchQuery(const CSeq_loc& sl, CScope& sc)

        CConstRef[CSeq_loc] GetQuerySeqLoc()
        CConstRef[CSeq_id] GetQueryId()
        CRef[CScope] GetScope()

        TSeqPos GetLength() const


    cppclass CBlastQueryVector(CObject):
        CBlastQueryVector()

        void AddQuery(CRef[CBlastSearchQuery] q)
        bool Empty()
        size_t Size()

        CRef[CBlastSearchQuery] operator[](size_t i)
