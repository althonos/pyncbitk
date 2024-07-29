from libcpp cimport bool

from ....corelib.ncbiobj cimport CObject, CRef
from ....objects.seqset.bioseq_set cimport CBioseq_set
from .blast_options cimport CBlastOptions


cdef extern from "algo/blast/api/query_data.hpp" namespace "ncbi::blast" nogil:

    cppclass ILocalQueryData(CObject):
        ILocalQueryData()

        # BLAST_SequenceBlk* GetSequenceBlk()
        # BlastQueryInfo* GetQueryInfo()
        size_t GetNumQueries()
        # CConstRef[CSeq_loc] GetSeq_loc(size_t index)
        size_t GetSeqLength(size_t index)
        size_t GetSumOfSequenceLengths()
        # void GetMessages(TSearchMessages& messages) const
        # void GetQueryMessages(size_t index, TQueryMessages& qmsgs)
        bool IsValidQuery(size_t index)
        bool IsAtLeastOneQueryValid()
        void FlushSequenceData()

    cppclass IRemoteQueryData(CObject):
        CRef[CBioseq_set] GetBioseqSet()
        # typedef list< CRef<objects::CSeq_loc> > TSeqLocs;
        # virtual TSeqLocs GetSeqLocs() = 0;

    cppclass IQueryFactory(CObject):
        CRef[ILocalQueryData] MakeLocalQueryData(const CBlastOptions* opts)
        CRef[IRemoteQueryData] MakeRemoteQueryData()