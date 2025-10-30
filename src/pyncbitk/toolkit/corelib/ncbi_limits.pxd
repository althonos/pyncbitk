from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from .ncbiype cimport Int1, Uint1, Int2, Uint2, Int4, Uint4, Int8, Uint8

cdef extern from "corelib/ncbi_limits.hpp" nogil:

    const Int1  kMin_I1
    const Int1  kMax_I1
    const Uint1 kMax_UI1

    const Int2  kMin_I2
    const Int2  kMax_I2
    const Uint2 kMax_UI2

    const Int4  kMin_I4
    const Int4  kMax_I4
    const Uint4 kMax_UI4

    const Int8  kMin_I8
    const Int8  kMax_I8
    const Uint8 kMax_UI8
    