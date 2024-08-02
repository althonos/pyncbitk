from libcpp cimport bool
from libcpp.vector cimport vector

from ...corelib.ncbimisc cimport TSeqPos
from ...serial.serialbase cimport CSerialObject, CStringAliasBase

cdef extern from "objects/seq/NCBI4na_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CNCBI4na_Base(CStringAliasBase[vector[char]]):
        CNCBI4na_Base()

cdef extern from "objects/seq/NCBI4na.hpp" namespace "ncbi::objects" nogil:

    cppclass CNCBI4na(CNCBI4na_Base):
        CNCBI4na()
        CNCBI4na(const vector[char]& data)
