from libcpp.vector cimport vector

from ....corelib.ncbiobj cimport CRef, CConstRef
from ....objects.seq.bioseq cimport CBioseq
from ....objects.seqset.bioseq_set cimport CBioseq_set
from .query_data cimport IQueryFactory
from .sseqloc cimport TSeqLocVector

cdef extern from "algo/blast/api/objmgrfree_query_data.hpp" namespace "ncbi::blast" nogil:

    cppclass CObjMgrFree_QueryFactory(IQueryFactory):
        CObjMgrFree_QueryFactory(CConstRef[CBioseq] bioseq)
        CObjMgrFree_QueryFactory(CConstRef[CBioseq_set] bioseq_set)
