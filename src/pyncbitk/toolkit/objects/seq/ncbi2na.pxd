from libcpp cimport bool
from libcpp.vector cimport vector

from ...serial.serialbase cimport CStringAliasBase

cdef extern from "objects/seq/NCBI2na_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CNCBI2na_Base(CStringAliasBase[vector[char]]):
        CNCBI2na_Base()

cdef extern from "objects/seq/NCBI2na.hpp" namespace "ncbi::objects" nogil:

    cppclass CNCBI2na(CNCBI2na_Base):
        CNCBI2na()
        CNCBI2na(const vector[char]& data)
