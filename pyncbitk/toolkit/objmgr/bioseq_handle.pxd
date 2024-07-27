# from libcpp cimport bool
# from libcpp.string cimport string
# from libcpp.vector cimport vector

# from ..corelib.ncbiobj cimport CObject, CRef
# from .object_manager cimport CObjectManager
from .scope cimport CScope

cdef extern from "objmgr/bioseq_handle.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_Handle:
        CBioseq_Handle()
        
        void Reset()
        CScope& GetScope() const

        