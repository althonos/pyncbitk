from ...serial.serialbase cimport CSerialObject

cdef extern from "objects/seqset/Bioseq_set_.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_set_Base(CSerialObject):
        CBioseq_set_Base()


cdef extern from "objects/seqset/Seq_entry.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_set(CBioseq_set_Base):
        CBioseq_set()