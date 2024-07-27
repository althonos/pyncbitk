from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t 

cdef extern from "corelib/ddumpable.hpp" namespace "ncbi" nogil:
    
    cppclass CDebugDumpable:
        CDebugDumpable()

