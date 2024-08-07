from libcpp cimport bool
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp.vector cimport vector

from ...corelib.ncbimisc cimport TSeqPos
from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from ..general.object_id cimport CObject_id
from ..seqloc.seq_loc cimport CSeq_loc
from ..seqloc.na_strand cimport ENa_strand
from ..seqloc.seq_id cimport CSeq_id
from .score cimport CScore


# cdef extern from "objects/seqalign/Seq_align_.hpp" namespace "ncbi::objects::CSeq_align_Base::C_Segs" nogil:

#     ctypedef CSerialObject Tparent



cdef extern from "objects/seqalign/Seq_align_.hpp" namespace "ncbi::objects::CSeq_align_Base" nogil:

    ctypedef CSerialObject Tparent
    
    enum EType:
        eType_not_set
        eType_global
        eType_diags
        eType_partial
        eType_disc
        eType_other

    cppclass C_Segs(CSerialObject):
        pass

    ctypedef EType TType
    ctypedef int TDim
    ctypedef vector[CRef[CScore]] TScore
    ctypedef C_Segs TSegs
    ctypedef list[CRef[CSeq_loc]] TBounds
    ctypedef list[CRef[CObject_id]] TId
    # ctypedef list[CRef[CUser_object]] TExt

    enum E_memberIndex:
        e__allMandatory
        e_type
        e_dim
        e_score
        e_segs
        e_bounds
        e_id
        e_ext

    

cdef extern from "objects/seqalign/Seq_align_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_align_Base(CSerialObject):
        CSeq_align_Base()

        bool IsSetType() const
        bool CanGetType() const
        void ResetType()
        TType GetType() const
        void SetType(TType value)
        TType& GetTypeMut "SetType" ()

        bool IsSetDim() const
        bool CanGetDim() const
        void ResetDim()
        TDim GetDim() const
        void SetDim(TDim value)
        TDim& GetDimMut "SetDim" ()

        bool IsSetScore() const
        bool CanGetScore() const
        void ResetScore()
        const TScore& GetScore() const
        TScore& GetScoreMut "SetScore" ()


cdef extern from "objects/seqalign/Seq_align.hpp" namespace "ncbi::objects::CSeq_align" nogil:

    ctypedef CSeq_align_Base Tparent
    ctypedef int TDim

    enum EScoreType:
        eScore_Score
        eScore_Blast
        eScore_BitScore
        eScore_EValue
        eScore_AlignLength
        eScore_IdentityCount
        eScore_PositiveCount
        eScore_NegativeCount
        eScore_MismatchCount
        eScore_GapCount
        eScore_PercentIdentity_Gapped
        eScore_PercentIdentity_Ungapped
        eScore_PercentIdentity_GapOpeningOnly
        eScore_PercentCoverage
        eScore_SumEValue
        eScore_CompAdjMethod
        eScore_HighQualityPercentCoverage
        eScore_Matches
        eScore_OverallIdentity
        eScore_Splices
        eScore_ConsensusSplices
        eScore_ProductCoverage
        eScore_ExonIdentity
        eScore_PercentIdentity

    struct SIndel:
        TSeqPos product_pos
        TSeqPos genomic_pos
        TDim row
        TSeqPos length


cdef extern from "objects/seqalign/Seq_align.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_align(CSeq_align_Base):
        CSeq_align()

        # CRange<TSeqPos> GetSeqRange(TDim row) const;
        TSeqPos         GetSeqStart(TDim row) except +
        TSeqPos         GetSeqStop (TDim row) except +
        ENa_strand      GetSeqStrand(TDim row) except +
        const CSeq_id&  GetSeq_id(TDim row) except +

        TSeqPos         GetTotalGapCount(TDim row = -1) except +
        # TSeqPos         GetTotalGapCountWithinRange(const TSeqRange &range, TDim row = -1) const
        # TSeqPos         GetTotalGapCountWithinRanges(const CRangeCollection<TSeqPos> &ranges, TDim row = -1) const

        TSeqPos         GetNumGapOpenings(TDim row = -1) except +
        # TSeqPos         GetNumGapOpeningsWithinRange(const TSeqRange &range, TDim row = -1) const;
        # TSeqPos         GetNumGapOpeningsWithinRanges(const CRangeCollection<TSeqPos> &ranges, TDim row = -1) const;

        TSeqPos         GetNumFrameshifts(TDim row = -1) except +
        # TSeqPos         GetNumFrameshiftsWithinRange(const TSeqRange &range, TDim row = -1) const;
        # TSeqPos         GetNumFrameshiftsWithinRanges(const CRangeCollection<TSeqPos> &ranges, TDim row = -1) const;

        vector[SIndel]  GetIndels(TDim row = -1) except +
        # vector[SIndel]  GetIndelsWithinRange(const TSeqRange &range, TDim row = -1) const;
        # vector[SIndel]  GetIndelsWithinRanges(const CRangeCollection<TSeqPos> &ranges, TDim row = -1) const;

        vector[SIndel]  GetFrameshifts(TDim row = -1) except +
        # vector[SIndel]  GetFrameshiftsWithinRange(const TSeqRange &range, TDim row = -1) const;
        # vector[SIndel]  GetFrameshiftsWithinRanges(const CRangeCollection<TSeqPos> &ranges,TDim row = -1) const;

        vector[SIndel]  GetNonFrameshifts(TDim row = -1) except +
        # vector[SIndel]  GetNonFrameshiftsWithinRange(const TSeqRange &range, TDim row = -1) const;
        # vector[SIndel]  GetNonFrameshiftsWithinRanges(const CRangeCollection<TSeqPos> &ranges, TDim row = -1) const;

        # CRangeCollection<TSeqPos> GetAlignedBases(TDim row) except +

        TSeqPos         GetAlignLength(bool include_gaps = true) except +

        # TSeqPos         GetAlignLengthWithinRange(const TSeqRange &range, bool include_gaps = true) const;
        # TSeqPos         GetAlignLengthWithinRanges(const CRangeCollection<TSeqPos> &ranges, bool include_gaps = true) const;

        double          AlignLengthRatio() const


        # bool GetNamedScore(const string& id, int &score) const
        # bool GetNamedScore(const string& id, double &score) const
        bool GetNamedScore(EScoreType type, int &score) const
        bool GetNamedScore(EScoreType type, double &score) const

        # void SetNamedScore(const string& id, int score)
        # void SetNamedScore(const string& id, double score)
        void SetNamedScore(EScoreType type, int score)
        void SetNamedScore(EScoreType type, double score)

        void ResetNamedScore(const string& name)
        void ResetNamedScore(EScoreType    type)