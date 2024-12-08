from libcpp.vector cimport vector

from ....corelib.ncbiobj cimport CRef
from ....objmgr.scope cimport CScope
from .query_data cimport IQueryFactory
from .sseqloc cimport TSeqLocVector, CBlastQueryVector

cdef extern from "algo/blast/api/objmgr_query_data.hpp" namespace "ncbi::blast" nogil:

    cppclass CObjMgr_QueryFactory(IQueryFactory):
        CObjMgr_QueryFactory(TSeqLocVector& queries)
        CObjMgr_QueryFactory(CBlastQueryVector& queries)

        vector[CRef[CScope]] ExtractScopes()
        # TSeqLocInfoVector ExtractUserSpecifiedMasks()
        TSeqLocVector GetTSeqLocVector()