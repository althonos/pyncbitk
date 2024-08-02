from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..corelib.ncbiobj cimport CObject, CRef
from .ncbi_tree cimport CTreePair
from .ncbistr cimport PEqualNocase_Conditional

cdef extern from "corelib/ncbi_config.hpp" namespace "ncbi::CConfig" nogil:

    ctypedef CTreePair[string, string, PEqualNocase_Conditional] TParamValue
    ctypedef TParamValue.TPairTreeNode TParamTree

cdef extern from "corelib/ncbi_config.hpp" namespace "ncbi" nogil:

    cppclass CConfig:
        pass