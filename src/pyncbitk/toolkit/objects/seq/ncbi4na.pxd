from libcpp cimport bool
from libcpp.vector cimport vector

from ...serial.serialbase cimport CStringAliasBase

cdef extern from "objects/seq/NCBI4na_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CNCBI4na_Base(CStringAliasBase[vector[char]]):
        CNCBI4na_Base()

cdef extern from "objects/seq/NCBI4na.hpp" namespace "ncbi::objects" nogil:

    cppclass CNCBI4na(CNCBI4na_Base):
        CNCBI4na()
        CNCBI4na(const vector[char]& data)
