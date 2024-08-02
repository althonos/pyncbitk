from libcpp.string cimport string

from .object_id cimport CObject_id
from ...serial.serialbase cimport CSerialObject

cdef extern from "objects/general/Dbtag_.hpp" namespace "ncbi::objects::CDbtag_Base" nogil:

    ctypedef string TDb
    ctypedef CObject_id TTag


cdef extern from "objects/general/Dbtag_.hpp" namespace "ncbi::objects" nogil:

    cppclass CDbtag_Base(CSerialObject):
        CDbtag_Base()

        bool IsSetDb(void) const
        bool CanGetDb(void) const
        void ResetDb(void)
        const TDb& GetDb(void) const
        void SetDb(const TDb& value)
        void SetDb(TDb&& value)
        TDb& GetDbMut "SetDb"(void)

        bool IsSetTag(void) const
        bool CanGetTag(void) const
        void ResetTag(void)
        const TTag& GetTag(void) const
        void SetTag(TTag& value)
        TTag& GetTagMut "SetTag"(void)

        virtual void Reset(void)

cdef extern from "objects/general/Dbtag.hpp" namespace "ncbi::objects" nogil:

    cppclass CDbtag(CDbtag_Base):
        CDbtag()

