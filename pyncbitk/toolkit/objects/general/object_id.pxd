from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CSerialObject

cdef extern from "objects/general/Object_id_.hpp" namespace "ncbi::objects::CObject_id_Base" nogil:

    enum E_Choice:
        e_not_set
        e_Id
        e_Str

    ctypedef int TId
    ctypedef string TStr


cdef extern from "objects/general/Object_id_.hpp" namespace "ncbi::objects" nogil:

    cppclass CObject_id_Base(CSerialObject):
        CObject_id_Base()

        void Reset()
        void ResetSelection()

        E_Choice Which()
        void CheckSelected(E_Choice index) const
        void ThrowInvalidSelection(E_Choice index) except +
        void Select(E_Choice index)
        # void Select(E_Choice index, EResetVariant reset)
        # void Select(E_Choice index, EResetVariant reset, CObjectMemoryPool* pool)

        bool IsId() const
        TId GetId() const
        TId& SetId()
        void SetId(TId value)

        # typedef string TStr
        bool IsStr() const
        const TStr& GetStr() const
        TStr& SetStr()
        void SetStr(const TStr& value)

cdef extern from "objects/general/Object_id.hpp" namespace "ncbi::objects" nogil:

    cppclass CObject_id(CObject_id_Base):
        CObject_id()

        bool Match(const CObject_id& oid2) const
        int Compare(const CObject_id& oid2) const
        