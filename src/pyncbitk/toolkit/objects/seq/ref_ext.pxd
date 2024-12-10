from ..seqloc.seq_loc cimport CSeq_loc


cdef extern from "objects/seq/Ref_ext_.hpp" namespace "ncbi::objects" nogil:

    cppclass CRef_ext_Base(CSeq_loc):
        CRef_ext_Base()

        const CSeq_loc& Get() const
        CSeq_loc& GetMut "Set"()


cdef extern from "objects/seq/Ref_ext.hpp" namespace "ncbi::objects" nogil:

    cppclass CRef_ext(CRef_ext_Base):
        CRef_ext()
