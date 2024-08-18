from libcpp cimport bool
from libcpp.vector cimport vector

from ...serial.serialbase cimport CStringAliasBase

cdef extern from "objects/seq/NCBI8na_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CNCBI8na_Base(CStringAliasBase[vector[char]]):
        CNCBI8na_Base()

cdef extern from "objects/seq/NCBI8na.hpp" namespace "ncbi::objects" nogil:

    cppclass CNCBI8na(CNCBI8na_Base):
        CNCBI8na()
        CNCBI8na(const vector[char]& data)
