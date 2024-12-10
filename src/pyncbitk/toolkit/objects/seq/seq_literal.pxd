from libcpp cimport bool

from ...serial.serialbase cimport CSerialObject
from ...corelib.ncbimisc cimport TSeqPos
from ..general.int_fuzz cimport CInt_fuzz
from .seq_data cimport CSeq_data

cdef extern from "objects/seq/Seq_literal_.hpp" namespace "ncbi::objects::CSeq_literal_Base" nogil:

    ctypedef TSeqPos TLength
    ctypedef CInt_fuzz TFuzz
    ctypedef CSeq_data TSeq_data

    enum E_memberIndex:
        e__allMandatory
        e_length
        e_fuzz
        e_seq_data


cdef extern from "objects/seq/Seq_literal_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_literal_Base(CSerialObject):
        CSeq_literal_Base()

        # mandatory
        bool IsSetLength() const
        bool CanGetLength() const
        void ResetLength()
        TLength GetLength() except +
        void SetLength(TLength value) except +
        TLength& GetLengthMut "SetLength"()

        # optional
        bool IsSetFuzz() const
        bool CanGetFuzz() const
        void ResetFuzz()
        const TFuzz& GetFuzz() except +
        void SetFuzz(TFuzz& value) except +
        TFuzz& GetFuzzMut "SetFuzz"()

        # optional
        bool IsSetSeq_data() const
        bool CanGetSeq_data() const
        void ResetSeq_data()
        const TSeq_data& GetSeq_data() except +
        void SetSeq_data(TSeq_data& value) except +
        TSeq_data& GetSeq_dataMut "SetSeq_data"()

        void Reset()



cdef extern from "objects/seq/Seq_literal.hpp" namespace "ncbi::objects::CSeq_literal" nogil:

    enum EBridgeableStatus:
        e_NotAGap
        e_Bridgeable
        e_NotBridgeable
        e_MaybeBridgeable


cdef extern from "objects/seq/Seq_literal.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_literal(CSeq_literal_Base):
        CSeq_literal()

        EBridgeableStatus GetBridgeability() const