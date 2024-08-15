from libcpp cimport bool

from ...serial.serialbase cimport CSerialObject, EResetVariant
from ...serial.serializable cimport CSerializable
from ..seq.bioseq cimport CBioseq
from .bioseq_set cimport CBioseq_set

cdef extern from "objects/seqset/Seq_entry_.hpp" namespace "ncbi::objects::CSeq_entry_Base" nogil:

    enum E_Choice:
        e_not_set
        e_Seq
        e_Set

    ctypedef CBioseq TSeq
    ctypedef CBioseq_set TSet


cdef extern from "objects/seqset/Seq_entry_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_entry_Base(CSerialObject):
        CSeq_entry_Base()

        void Reset()
        void ResetSelection()

        E_Choice Which() const
        void CheckSelected(E_Choice index) const

        void Select(E_Choice index)
        void Select(E_Choice index, EResetVariant reset)
        # void Select(E_Choice index, EResetVariant reset, CObjectMemoryPool* pool)

        bool IsSeq() const
        const TSeq& GetSeq() const
        TSeq& GetSeqMut "SetSeq"()
        void SetSeq(TSeq& value)

        bool IsSet() const
        const TSet& GetSet() const
        TSet& GetSetMut "SetSet"()
        void SetSet(TSet& value)


cdef extern from "objects/seqset/Seq_entry.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_entry(CSeq_entry_Base, CSerializable):
        CSeq_id()