from libcpp cimport bool

from ...serial.serialbase cimport CSerialObject
from ...objects.general.object_id cimport CObject_id

cdef extern from "objects/seqalign/Score_.hpp" namespace "ncbi::objects::CScore_Base::C_Value" nogil:

    enum E_Choice:
        e_not_set
        e_Real
        e_Int

    ctypedef double TReal
    ctypedef int TInt


cdef extern from "objects/seqalign/Score_.hpp" namespace "ncbi::objects::CScore_Base" nogil:

    cppclass C_Value(CSerialObject):
        C_Value()

        void Reset()
        void ResetSelection()
        E_Choice Which() const
        void CheckSelected(E_Choice index) except +

        bool IsReal() const
        TReal GetReal() const
        TReal& SetReal()
        void SetReal(TReal value)
    
        bool IsInt() const
        TInt GetInt() const
        TInt& SetInt()
        void SetInt(TInt value)

    ctypedef C_Value TValue
    ctypedef CObject_id TId


cdef extern from "objects/seqalign/Score_.hpp" namespace "ncbi::objects" nogil:

    cppclass CScore_Base(CSerialObject):
        CScore_Base()

        bool IsSetId() const
        bool CanGetId() const
        void ResetId()
        const TId& GetId() const
        void SetId(TId& value)
        TId& GetIdMut "SetId" ()

        bool IsSetValue() const
        bool CanGetValue() const
        void ResetValue()
        const TValue& GetValue() const
        TValue& GetValueMut "SetValue"()

        void Reset()
        void DoNotDeleteThisObject()

cdef extern from "objects/seqalign/Score.hpp" namespace "ncbi::objects" nogil:

    cppclass CScore(CScore_Base):
        CScore()
