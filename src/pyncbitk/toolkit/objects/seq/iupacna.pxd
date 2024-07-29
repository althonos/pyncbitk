from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.ncbimisc cimport TSeqPos
from ...serial.serialbase cimport CSerialObject, CStringAliasBase

cdef extern from "objects/seq/IUPACna_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CIUPACna_Base(CStringAliasBase[string]):
        CIUPACna_Base()

cdef extern from "objects/seq/IUPACna.hpp" namespace "ncbi::objects" nogil:

    cppclass CIUPACna(CIUPACna_Base):
        CIUPACna()
        CIUPACna(string& data)
