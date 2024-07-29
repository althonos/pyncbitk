from libcpp cimport bool
from libcpp.list cimport list

from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from .seq_align cimport CSeq_align


cdef extern from "objects/seqalign/Seq_align_set_.hpp" namespace "ncbi::objects::CSeq_align_set_Base" nogil:

    # ctypedef CSerialObject Tparent
    ctypedef list[CRef[CSeq_align]] Tdata
        

cdef extern from "objects/seqalign/Seq_align_set_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_align_set_Base(CSerialObject):
        CSeq_align_set_Base()

        bool IsSet() const
        bool CanGet() const
        void Reset()
        const Tdata& Get() const
        Tdata& GetMut "Set" ()


cdef extern from "objects/seqalign/Seq_align_set.hpp" namespace "ncbi::objects::CSeq_align_set" nogil:

    # ctypedef CSeq_align_set_Base Tparent
    ctypedef int TDim


cdef extern from "objects/seqalign/Seq_align_set.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_align_set(CSeq_align_set_Base):
        CSeq_align_set()

        void SwapRows(TDim row1, TDim row2)
        Tdata.size_type Size() const
        bool IsEmpty() const