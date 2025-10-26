from libcpp cimport bool
from libcpp.string cimport string
from libcpp.list cimport list

from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from .seqdesc cimport CSeqdesc


cdef extern from "objects/seq/Seq_descr_.hpp" namespace "ncbi::objects::CSeq_descr_Base" nogil:

    ctypedef list[CRef[CSeqdesc]] Tdata


cdef extern from "objects/seq/Seq_descr_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_descr_Base(CSerialObject):
        CSeq_descr_Base()

        bool IsSet() const
        bool CanGet() const
        void Reset()
        const Tdata& Get() except +
        Tdata& GetMut "Set"() except +

     
cdef extern from "objects/seq/Seq_descr.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_descr(CSeq_descr_Base):
        CSeq_descr()
    