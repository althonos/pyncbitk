from libcpp cimport bool
from libcpp.vector cimport vector

from ....corelib.ncbitype cimport Uint4 
from ....corelib.ncbiobj cimport CRef, CConstRef 
from ....objects.seqloc.seq_loc cimport CSeq_loc
from ....objmgr.scope cimport CScope

cdef extern from "algo/blast/api/sseqloc.hpp" namespace "ncbi::blast" nogil:

    cppclass SSeqLoc:
        CConstRef[CSeq_loc] seqloc
        CRef[CScope] score
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