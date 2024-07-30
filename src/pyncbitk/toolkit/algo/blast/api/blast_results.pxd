from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

from ....corelib.ncbitype cimport Int8
from ....corelib.ncbiobj cimport CObject, CRef, CConstRef
from ....objects.seqloc.seq_id cimport CSeq_id
from ....objects.seqalign.seq_align_set cimport CSeq_align_set
from .setup_factory cimport CThreadable
from .blast_types cimport EResultType, TSeqAlignVector, TSearchMessages



cdef extern from "algo/blast/api/blast_results.hpp" namespace "ncbi::blast" nogil:

    cppclass CBlastAncillaryData(CObject):
        # CBlastAncillaryData(EBlastProgramType, int query_number, const BlastScoreBlk* sbp, const BlastQueryInfo* query_info)
        CBlastAncillaryData(pair[double, double] lambda_, pair[double, double] k, pair[double, double] h, Int8 effective_search_space)
        CBlastAncillaryData(pair[double, double] lambda_, pair[double, double] k, pair[double, double] h, Int8 effective_search_space, bool is_psiblast)

        CBlastAncillaryData(const CBlastAncillaryData& rhs)
        CBlastAncillaryData& operator=(const CBlastAncillaryData& rhs)

        # const Blast_GumbelBlk * GetGumbelBlk() const
        # const Blast_KarlinBlk * GetUngappedKarlinBlk() const
        # const Blast_KarlinBlk * GetGappedKarlinBlk() const
        # const Blast_KarlinBlk * GetPsiUngappedKarlinBlk() const 
        # const Blast_KarlinBlk * GetPsiGappedKarlinBlk() const
        Int8 GetSearchSpace() const
        Int8 GetLengthAdjustment() const
        void SetLengthAdjustment(int len_adj)

    cppclass CSearchResults(CObject):
        # CSearchResults(
        #     CConstRef[CSeq_id]          query,
        #     CRef[CSeq_align_set]        align, 
        #     const TQueryMessages      & errs,
        #     CRef[CBlastAncillaryData]   ancillary_data,
        #     const TMaskedQueryRegions * query_masks = NULL,
        #     const string              & rid = kEmptyStr,
        #     const SPHIQueryInfo       * phi_query_info = NULL
        # )

        string GetRID() const
        void SetRID(const string& rid) noexcept

        CConstRef[CSeq_align_set] GetSeqAlign() const
        CRef[CSeq_align_set] GetSeqAlignMut "SetSeqAlign" () const

        bool HasAlignments() const

        CConstRef[CSeq_id] GetSeqId() const
        CRef[CBlastAncillaryData] GetAncillaryData() const

        # TQueryMessages GetErrors()
        # TQueryMessages GetErrors(int min_severity = eBlastSevError)

        bool HasErrors() const
        bool HasWarnings() const

        string GetErrorStrings() const
        string GetWarningStrings() const

        # void GetMaskedQueryRegions(TMaskedQueryRegions& flt_query_regions) const
        # void SetMaskedQueryRegions(const TMaskedQueryRegions& flt_query_regions)

        # void GetSubjectMasks(TSeqLocInfoVector& subj_masks) const
        # void SetSubjectMasks(const TSeqLocInfoVector& subj_masks)

        # const SPHIQueryInfo * GetPhiQueryInfo() const

        # void TrimSeqAlign(objects::CSeq_align_set::Tdata::size_type max_size)

cdef extern from "algo/blast/api/blast_results.hpp" namespace "ncbi::blast::CSearchResultSet" nogil:

    ctypedef CRef[CSearchResults] value_type
    ctypedef vector[CConstRef[CSeq_id]] TQueryIdVector
    ctypedef vector[value_type].size_type size_type
    ctypedef vector[CRef[CBlastAncillaryData]] TAncillaryVector
    ctypedef vector[value_type].const_iterator const_iterator
    ctypedef vector[value_type].iterator iterator

cdef extern from "algo/blast/api/blast_results.hpp" namespace "ncbi::blast" nogil:

    cppclass CSearchResultSet(CObject):
        CSearchResultSet()
        CSearchResultSet(EResultType res_type)

        CSearchResultSet(TSeqAlignVector aligns, TSearchMessages msg_vec)
        CSearchResultSet(TSeqAlignVector aligns, TSearchMessages msg_vec, EResultType res_type)

        CSearchResults& operator[](size_type i)
        # const CSearchResults& operator[](size_type i) const

        CSearchResults& GetResults(size_type qi, size_type si)
        # const CSearchResults& GetResults(size_type qi, size_type si) const

        CRef[CSearchResults] operator[](const CSeq_id & ident)
        # CConstRef[CSearchResults] operator[](const CSeq_id & ident) const

        size_type GetNumResults() const
        size_type GetNumQueries() const

        # void SetFilteredQueryRegions(const TSeqLocInfoVector& masks)
        # TSeqLocInfoVector GetFilteredQueryRegions() const

        size_type size() const 
        bool empty() const 
        const_iterator begin() const
        const_iterator end() const
        iterator begin() except +
        iterator end() except +
        void clear() noexcept
        void push_back(value_type& element) except +
        EResultType GetResultType() const
        void SetRID(const string& rid) except +