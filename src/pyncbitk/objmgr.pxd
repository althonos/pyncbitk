# cython: language_level=3, linetrace=True, binding=True

from .toolkit.corelib.ncbiobj cimport CRef
from .toolkit.objmgr.object_manager cimport CObjectManager
from .toolkit.objmgr.bioseq_handle cimport CBioseq_Handle
from .toolkit.objmgr.scope cimport CScope

from .objects.seq cimport BioSeq


cdef class BioSeqHandle:
    cdef CBioseq_Handle _handle


cdef class ObjectManager:
    cdef CRef[CObjectManager] _mgr

    cpdef Scope scope(self)


cdef class Scope:
    cdef CRef[CScope] _scope

    cpdef BioSeqHandle add_bioseq(self, BioSeq seq) except *