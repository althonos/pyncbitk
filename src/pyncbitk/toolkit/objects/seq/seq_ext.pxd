from libcpp cimport bool
from libcpp.string cimport string

from ..seqloc.seq_loc cimport CSeq_loc
from ...serial.serialbase cimport CSerialObject
from .ref_ext cimport CRef_ext
from .delta_ext cimport CDelta_ext


cdef extern from "objects/seq/Seq_ext_.hpp" namespace "ncbi::objects::CSeq_ext_Base" nogil:

    enum E_Choice:
        e_not_set
        e_Seg
        e_Ref
        e_Map
        e_Delta

    enum E_ChoiceStopper:
        e_MaxChoice

    # ctypedef CSeg_ext TSeg
    ctypedef CRef_ext TRef
    # ctypedef CMap_ext TMap
    ctypedef CDelta_ext TDelta


cdef extern from "objects/seq/Seq_ext_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_ext_Base(CSerialObject):
        CSeq_ext_Base()

        void Reset()
        void ResetSelection()

        E_Choice Which() const
        void CheckSelected(E_Choice index) const
        void ThrowInvalidSelection(E_Choice index) except +

        @staticmethod
        string SelectionName(E_Choice index)

        void Select(E_Choice index)
        # void Select(E_Choice index, EResetVariant reset)
        # void Select(E_Choice index, EResetVariant reset, CObjectMemoryPool* pool)

        # bool IsSeg() const
        # const TSeg& GetSeg() const
        # TSeg& GetSegMut "SetSeg"()
        # void SetSeg(TSeg& value)

        bool IsRef() const
        const TRef& GetRef() except +
        TRef& GetRefMut "SetRef"() except +
        void SetRef(TRef& value) except +

        # bool IsMap() const
        # const TMap& GetMap() const
        # TMap& SetMap()
        # void SetMap(TMap& value)

        bool IsDelta() const
        const TDelta& GetDelta() except +
        TDelta& GetDeltaMut "SetDelta"() except +
        void SetDelta(TDelta& value) except +
     
cdef extern from "objects/seq/Seq_ext.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_ext(CSeq_ext_Base):
        CSeq_ext()
    