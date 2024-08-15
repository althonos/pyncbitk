from libcpp cimport bool

from ...serial.serialbase cimport CSerialObject
from .seq_id cimport CSeq_id
from .na_strand cimport ENa_strand
from .seq_interval cimport CSeq_interval

cdef extern from "objects/seqloc/Seq_loc_.hpp" namespace "ncbi::objects::CSeq_loc_Base" nogil:

    enum E_Choice:
        e_not_set
        e_Null
        e_Empty
        e_Whole
        e_Int
        e_Packed_int
        e_Pnt
        e_Packed_pnt
        e_Mix
        e_Equiv
        e_Bond
        e_Feat

    ctypedef CSeq_id TEmpty
    ctypedef CSeq_id TWhole
    ctypedef CSeq_interval TInt
    # typedef CPacked_seqint TPacked_int
    # typedef CSeq_point TPnt
    # typedef CPacked_seqpnt TPacked_pnt
    # typedef CSeq_loc_mix TMix
    # typedef CSeq_loc_equiv TEquiv
    # typedef CSeq_bond TBond
    # typedef CFeat_id TFeat


cdef extern from "objects/seqloc/Seq_loc_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_loc_Base(CSerialObject):
        CSeq_loc_Base()

        void Reset()
        void ResetSelection()

        E_Choice Which() const
        void CheckSelected(E_Choice index) except + 
        void ThrowInvalidSelection(E_Choice index) except + 

        void Select(E_Choice index)
        # void Select(E_Choice index, EResetVariant reset)
        # void Select(E_Choice index, EResetVariant reset, CObjectMemoryPool* pool)

        bool IsNull() const
        void SetNull()

        bool IsEmpty() const
        const TEmpty& GetEmpty() const
        TEmpty& SetEmpty()
        void SetEmpty(TEmpty& value)

        bool IsWhole() const
        const TWhole& GetWhole() const
        TWhole& GetWholeMut "SetWhole" ()
        void SetWhole(TWhole& value)

        bool IsInt() const
        const TInt& GetInt() const
        TInt& GetIntMut "SetInt"()
        void SetInt(TInt& value)

        # bool IsPacked_int() const
        # const TPacked_int& GetPacked_int() const
        # TPacked_int& SetPacked_int()
        # void SetPacked_int(TPacked_int& value)

        # bool IsPnt() const
        # const TPnt& GetPnt() const
        # TPnt& SetPnt()
        # void SetPnt(TPnt& value)

        # bool IsPacked_pnt() const
        # const TPacked_pnt& GetPacked_pnt() const
        # TPacked_pnt& SetPacked_pnt()
        # void SetPacked_pnt(TPacked_pnt& value)

        # bool IsMix() const
        # const TMix& GetMix() const
        # TMix& SetMix()
        # void SetMix(TMix& value)

        # bool IsEquiv() const
        # const TEquiv& GetEquiv() const
        # TEquiv& SetEquiv()
        # void SetEquiv(TEquiv& value)

        # bool IsBond() const
        # const TBond& GetBond() const
        # TBond& SetBond()
        # void SetBond(TBond& value)

        # bool IsFeat() const
        # const TFeat& GetFeat() const
        # TFeat& SetFeat()
        # void SetFeat(TFeat& value)


cdef extern from "objects/seqloc/Seq_loc.hpp" namespace "ncbi::objects::CSeq_loc" nogil:
    pass
    # typedef CSeq_loc_Base Tparent;
    # typedef CPacked_seqpnt_Base::TPoints TPoints;
    # typedef CPacked_seqint_Base::Tdata   TIntervals;
    # typedef CSeq_loc_mix_Base::Tdata     TLocations;
    # typedef CSeq_id                      TId;
    # typedef ENa_strand                   TStrand;
    # typedef TSeqPos                      TPoint;
    # typedef CPacked_seqint::TRanges      TRanges;

cdef extern from "objects/seqloc/Seq_loc.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_loc(CSeq_loc_Base):
        CSeq_loc()
        CSeq_loc(E_Choice index)
        # CSeq_loc(TId& id, TPoint point)
        # CSeq_loc(TId& id, TPoint point, TStrand strand)
        # CSeq_loc(TId& id, const TPoints& points)
        # CSeq_loc(TId& id, const TPoints& points, TStrand strand)
        # CSeq_loc(TId& id, TPoint from_, TPoint to)
        # CSeq_loc(TId& id, TPoint from_, TPoint to, TStrand strand)
        # CSeq_loc(TId& id, TRanges ivals)
        # CSeq_loc(TId& id, TRanges ivals, TStrand strand)

        # bool IsSetStrand() const
        # bool IsSetStrand(EIsSetStrand flag) const
        # ENa_strand GetStrand() const
        # bool IsReverseStrand() const
        # void FlipStrand()
        # void SetStrand(ENa_strand strand)
        # void ResetStrand()

        # TSeqPos GetStart(ESeqLocExtremes ext) const
        # TSeqPos GetStop (ESeqLocExtremes ext) const

        # TSeqPos GetCircularLength(TSeqPos seq_len) const

        # void GetLabel(string* label) const

        # bool IsPartialStart(ESeqLocExtremes ext) const
        # bool IsPartialStop(ESeqLocExtremes ext) const
        
        # void SetPartialStart(bool val, ESeqLocExtremes ext)
        # void SetPartialStop (bool val, ESeqLocExtremes ext)

        # bool IsTruncatedStart(ESeqLocExtremes ext) const
        # bool IsTruncatedStop (ESeqLocExtremes ext) const

        # void SetTruncatedStart(bool val, ESeqLocExtremes ext)
        # void SetTruncatedStop (bool val, ESeqLocExtremes ext)

        # const CSeq_id* GetId(void) const

        # bool CheckId(const CSeq_id*& id, bool may_throw = true) const
        # void InvalidateIdCache(void) const

        # void SetId(CSeq_id& id)
        # void SetId(const CSeq_id& id)

        void InvalidateCache() const

        