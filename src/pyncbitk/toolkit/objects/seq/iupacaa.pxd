from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CStringAliasBase

cdef extern from "objects/seq/IUPACaa_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CIUPACaa_Base(CStringAliasBase[string]):
        CIUPACaa_Base()

cdef extern from "objects/seq/IUPACaa.hpp" namespace "ncbi::objects" nogil:

    cppclass CIUPACaa(CIUPACaa_Base):
        CIUPACaa()
        CIUPACaa(string& data)
