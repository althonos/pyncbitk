from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CSerialObject

cdef extern from "objects/general/Int_fuzz_.hpp" namespace "ncbi::objects::CInt_fuzz_Base" nogil:

    pass


cdef extern from "objects/general/Int_fuzz_.hpp" namespace "ncbi::objects" nogil:

    cppclass CInt_fuzz_Base(CSerialObject):
        CInt_fuzz_Base()



cdef extern from "objects/general/Int_fuzz.hpp" namespace "ncbi::objects" nogil:

    cppclass CInt_fuzz(CInt_fuzz_Base):
        CInt_fuzz()

      