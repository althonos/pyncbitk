from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..corelib.ncbiobj cimport CObject, CRef
from .ncbi_config cimport TParamTree

cdef extern from "corelib/plugin_manager.hpp" namespace "ncbi" nogil:

    ctypedef TParamTree TPluginManagerParamTree

    cppclass CPluginManagerBase(CObject):
        pass

    cppclass CPluginManager(CPluginManagerBase):
        pass