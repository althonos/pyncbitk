from libcpp.string cimport string

cdef extern from "corelib/ncbistr.hpp" nogil:

    const char* const kEmptyCStr
    const string& kEmptyStr