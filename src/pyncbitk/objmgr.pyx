# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from .toolkit.objmgr.object_manager cimport CObjectManager
from .toolkit.objmgr.scope cimport CScope

from .objects.seq cimport BioSeq


cdef class ObjectManager:
    def __init__(self):
        self._mgr = CObjectManager.GetInstance()

    cpdef Scope scope(self):
        return Scope(self)

    def register_data_loader(self, str name not None):
        cdef bytes _name = name.encode()
        self._mgr.GetObject().RegisterDataLoader(NULL, <string> _name)

    def get_registered_names(self):
        cdef vector[string] names
        self._mgr.GetObject().GetRegisteredNames(names)
        return list(names)


cdef class Scope:
    def __init__(self, ObjectManager manager):
        self._scope.Reset(new CScope(manager._mgr.GetObject()))

    def __enter__(self):
        return self

    def __exit__(self, exc_ty, exc_val, traceback):
        self.close()
        return False

    def close(self):
        self._scope.ReleaseOrNull()

    cpdef void add_bioseq(self, BioSeq seq) except *:
        self._scope.GetObject().AddBioseq(seq._ref.GetObject())