from libcpp.string cimport string

cdef extern from "corelib/ncbistr.hpp" nogil:

    const char* const kEmptyCStr
    const string& kEmptyStr


    cppclass PNocase_Conditional_Generic[T]:
        pass

    ctypedef PNocase_Conditional_Generic[string]      PNocase_Conditional
    ctypedef PNocase_Conditional_Generic[const char*] PNocase_Conditional_CStr

    cppclass PEqualNocase_Conditional_Generic[T](PNocase_Conditional_Generic[T]):
        pass

    ctypedef PEqualNocase_Conditional_Generic[string]       PEqualNocase_Conditional
    ctypedef PEqualNocase_Conditional_Generic[const char *] PEqualNocase_Conditional_CStr
