from libcpp cimport bool
from libcpp.string cimport string

from ....corelib.ncbiobj cimport CObject, CRef, CConstRef
from .setup_factory cimport CThreadable
from .blast_results cimport CSearchResultSet
from .blast_options_handle cimport CBlastOptionsHandle
from .uniform_search cimport CSearchDatabase
from .query_data cimport IQueryFactory


cdef extern from "algo/blast/api/local_db_adapter.hpp" namespace "ncbi::blast" nogil:

    cppclass CLocalDbAdapter(CObject):
        CLocalDbAdapter(const CSearchDatabase& dbinfo)
        CLocalDbAdapter(CRef[IQueryFactory] subject_sequences, CRef[CBlastOptionsHandle] opts_handle) except +
        CLocalDbAdapter(CRef[IQueryFactory] subject_sequences, CRef[CBlastOptionsHandle] opts_handle, bool dbscan_mode) except +
        CLocalDbAdapter(CRef[IQueryFactory] subject_sequences, CConstRef[CBlastOptionsHandle] opts_handle) except +
        CLocalDbAdapter(CRef[IQueryFactory] subject_sequences, CConstRef[CBlastOptionsHandle] opts_handle, bool dbscan_mode) except +
        # CLocalDbAdapter(BlastSeqSrc* seqSrc, CRef[IBlastSeqInfoSrc] seqInfoSrc)

        void ResetBlastSeqSrcIteration()
        # BlastSeqSrc* MakeSeqSrc()
        # IBlastSeqInfoSrc* MakeSeqInfoSrc()
        int GetFilteringAlgorithm()
        string GetFilteringAlgorithmKey()

        bool IsProtein() const
        string GetDatabaseName() const
        bool IsBlastDb() const
        bool IsDbScanMode() const
        CRef[CSearchDatabase] GetSearchDatabase() const