from libcpp cimport bool

from ...corelib.ncbimisc cimport TSeqPos
from ...serial.serialbase cimport CSerialObject
from ..general.int_fuzz cimport CInt_fuzz
from .na_strand cimport ENa_strand, ESeqLocExtremes
from .seq_id cimport CSeq_id

cdef extern from "objects/seqloc/Seq_interval_.hpp" namespace "ncbi::objects::Seq_interval_" nogil:
    ctypedef TSeqPos TFrom
    ctypedef TSeqPos TTo
    ctypedef ENa_strand TStrand
    ctypedef CSeq_id TId
    ctypedef CInt_fuzz TFuzz_from
    ctypedef CInt_fuzz TFuzz_to

cdef extern from "objects/seqloc/Seq_interval_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_interval_Base(CSerialObject):
        CSeq_interval_Base()

        bool IsSetFrom() const
        bool CanGetFrom() const
        void ResetFrom()
        TFrom GetFrom() const
        void SetFrom(TFrom value)
        TFrom& GetFromMut "SetFrom"()

        bool IsSetTo() const
        bool CanGetTo() const
        void ResetTo()
        TTo GetTo() const
        void SetTo(TTo value);
        TTo& GetToMut "SetTo"()

        bool IsSetStrand() const
        bool CanGetStrand() const
        void ResetStrand()
        TStrand GetStrand() const
        void SetStrand(TStrand value)
        TStrand& GetStrandMut "SetStrand"()

        bool IsSetId() const
        bool CanGetId() const
        void ResetId()
        const TId& GetId() const
        void SetId(TId& value)
        TId& GetIdMut "SetId"()

        bool IsSetFuzz_from() const
        bool CanGetFuzz_from() const
        void ResetFuzz_from()
        const TFuzz_from& GetFuzz_from() const
        void SetFuzz_from(TFuzz_from& value)
        TFuzz_from& GetFuzz_fromMut "SetFuzz_from"()

        bool IsSetFuzz_to() const
        bool CanGetFuzz_to() const
        void ResetFuzz_to()
        const TFuzz_to& GetFuzz_to() const
        void SetFuzz_to(TFuzz_to& value)
        TFuzz_to& GetFuzz_toMut "SetFuzz_to"()

        void Reset()

cdef extern from "objects/seqloc/Seq_interval.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_interval(CSeq_interval_Base):
        CSeq_interval()
        CSeq_interval(TId& id, TFrom from_, TTo to)
        CSeq_interval(TId& id, TFrom from_, TTo to, TStrand strand)

        TSeqPos GetLength() const

        bool IsPartialStart(ESeqLocExtremes ext) const
        bool IsPartialStop (ESeqLocExtremes ext) const

        void SetPartialStart(bool val, ESeqLocExtremes ext)
        void SetPartialStop (bool val, ESeqLocExtremes ext)

        bool IsTruncatedStart(ESeqLocExtremes ext) const
        bool IsTruncatedStop (ESeqLocExtremes ext) const

        void SetTruncatedStart(bool val, ESeqLocExtremes ext)
        void SetTruncatedStop (bool val, ESeqLocExtremes ext)

        TSeqPos GetStart(ESeqLocExtremes ext) const
        TSeqPos GetStop (ESeqLocExtremes ext) const

        void FlipStrand()