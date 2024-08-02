from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CStringAliasBase

cdef extern from "objects/seq/NCBIeaa_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CNCBIeaa_Base(CStringAliasBase[string]):
        CNCBIeaa_Base()

cdef extern from "objects/seq/NCBI4na.hpp" namespace "ncbi::objects" nogil:

    cppclass CNCBIeaa(CNCBIeaa_Base):
        CNCBIeaa()
        CNCBIeaa(const string& data)
