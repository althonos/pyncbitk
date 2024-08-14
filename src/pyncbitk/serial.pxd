# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.string cimport string

from .toolkit.serial.serialbase cimport CSerialObject

cdef class Serial:
    
    cdef CSerialObject* _serial(self)

    cpdef string dumps(self, str format=*, bool indent=*, bool eol=*)