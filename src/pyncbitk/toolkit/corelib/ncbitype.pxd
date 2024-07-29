from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t 

cdef extern from "corelib/ncbitype.h" nogil:

    ctypedef char Char
    ctypedef signed char Schar
    ctypedef unsigned char Uchar

    ctypedef int8_t Int1
    ctypedef uint8_t Uint1
    ctypedef int16_t Int2
    ctypedef uint16_t Uint2
    ctypedef int32_t Int4
    ctypedef uint32_t Uint4
    ctypedef int64_t Int8
    ctypedef uint64_t Uint8

