from ...serial.serialbase cimport CSerialObject


cdef extern from "objects/seqalign/Dense_diag_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CDense_diag_Base(CSerialObject):
        CDense_diag_Base()


cdef extern from "objects/seqalign/Dense_diag.hpp" namespace "ncbi::objects" nogil:

    cppclass CDense_diag(CDense_diag_Base):
        CDense_diag()