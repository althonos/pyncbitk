from ....corelib.ncbiobj cimport CObject, CRef
from .setup_factory cimport CThreadable
from .blast_results cimport CSearchResultSet
from .blast_types cimport TSearchMessages
from .blast_options_handle cimport CBlastOptionsHandle
from .uniform_search cimport CSearchDatabase
from .query_data cimport IQueryFactory
from .local_db_adapter cimport CLocalDbAdapter


cdef extern from "algo/blast/api/local_blast.hpp" namespace "ncbi::blast" nogil:

    cppclass CLocalBlast(CObject, CThreadable):
        CLocalBlast(CRef[IQueryFactory] query_factory, CRef[CBlastOptionsHandle] opts_handle, const CSearchDatabase& dbinfo) except +
        CLocalBlast(CRef[IQueryFactory] query_factory, CRef[CBlastOptionsHandle] opts_handle, CRef[CLocalDbAdapter] db) except +
        # CLocalBlast(CRef[IQueryFactory] query_factory, CRef[CBlastOptionsHandle] opts_handle, BlastSeqSrc* seqsrc, CRef[IBlastSeqInfoSrc] seqInfoSrc)

        CRef[CSearchResultSet] Run() except +

        TSearchMessages GetSearchMessages() const
        # Int4 GetNumExtensions()
        # BlastDiagnostics* GetDiagnostics()
        # void SetBatchNumber( int batch_num )