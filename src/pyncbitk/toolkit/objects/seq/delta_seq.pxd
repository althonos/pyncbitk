from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CSerialObject
from ..seqloc.seq_loc cimport CSeq_loc
from .seq_literal cimport CSeq_literal


cdef extern from "objects/seq/Delta_seq_.hpp" namespace "ncbi::objects::CDelta_seq_Base":

    enum E_Choice:
        e_not_set
        e_Loc
        e_Literal

    enum E_ChoiceStopper:
        e_MaxChoice

    ctypedef CSeq_loc TLoc
    ctypedef CSeq_literal TLiteral



cdef extern from "objects/seq/Delta_seq_.hpp" namespace "ncbi::objects":

    cppclass CDelta_seq_Base(CSerialObject):
        CDelta_seq_Base()

        void Reset()
        void ResetSelection()

        E_Choice Which() const

        void CheckSelected(E_Choice index) except +
        void ThrowInvalidSelection(E_Choice index) except +

        @staticmethod
        string SelectionName(E_Choice index) except +

        void Select(E_Choice index) except +
        # void Select(E_Choice index, EResetVariant reset) except +
        # void Select(E_Choice index, EResetVariant reset, CObjectMemoryPool* pool) except +

        bool IsLoc() const
        const TLoc& GetLoc() except +
        TLoc& GetLocMut "SetLoc"() except +
        void SetLoc(TLoc& value)

        bool IsLiteral() const
        const TLiteral& GetLiteral() except +
        TLiteral& GetLiteralMut "SetLiteral"() except +
        void SetLiteral(TLiteral& value) except +


cdef extern from "objects/seq/Delta_seq.hpp" namespace "ncbi::objects":

    cppclass CDelta_seq(CDelta_seq_Base):
        CDelta_seq()
