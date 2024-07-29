from libcpp.string cimport string

cdef extern from "corelib/tempstr.hpp" namespace "ncbi::CTempString" nogil:

    ctypedef char        value_type
    ctypedef size_t      size_type
    ctypedef const char* const_iterator


cdef extern from "corelib/tempstr.hpp" namespace "ncbi" nogil:

    cppclass CTempString:
        CTempString()
        CTempString(const char* str)
        CTempString(const char* str, size_type len)

        CTempString(const string& str)
        CTempString(const string& str, size_type pos, size_type len)

        CTempString(const CTempString& str)
        CTempString(const CTempString& str, size_type pos)
        CTempString(const CTempString& str, size_type pos, size_type len)
