from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..corelib.ncbiobj cimport CObject, CRef

cdef extern from "objmgr/object_manager.hpp" namespace "ncbi::objects::CObjectManager" nogil:

    enum EIsDefault:
        eDefault
        eNonDefault

    ctypedef int TPriority
    enum EPriority:
        kPriority_Entry
        kPriority_Local
        kPriority_Replace
        kPriority_Loader
        kPriority_Extra
        kPriority_Default
        kPriority_NotSet

    ctypedef vector[string] TRegisteredNames




cdef extern from "objmgr/object_manager.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CObjectManager(CObject):
        @staticmethod
        CRef[CObjectManager] GetInstance()

        # CDataLoader* RegisterDataLoader(TPluginManagerParamTree* params = 0, const string& driver_name = kEmptyStr);
        # CDataLoader* FindDataLoader(const string& loader_name) const
        void GetRegisteredNames(TRegisteredNames& names)
        # void SetLoaderOptions(const string& loader_name,
        #                   EIsDefault    is_default,
        #                   TPriority     priority = kPriority_Default);
        # bool RevokeDataLoader(CDataLoader& loader)
        bool RevokeDataLoader(const string& loader_name)
        void RevokeAllDataLoaders()