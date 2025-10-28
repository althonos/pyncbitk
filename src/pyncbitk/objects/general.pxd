# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.general.dbtag cimport CDbtag
from ..toolkit.objects.general.object_id cimport CObject_id

from ..serial cimport Serial

# --- ObjectId -----------------------------------------------------------------

cdef class ObjectId(Serial):
    cdef CRef[CObject_id] _ref

    @staticmethod
    cdef ObjectId _wrap(CRef[CObject_id] ref)

# --- DBTag --------------------------------------------------------------------

cdef class DBTag(Serial):
    cdef CRef[CDbtag] _ref

    @staticmethod
    cdef DBTag _wrap(CRef[CDbtag] ref)

    cpdef str url(self, object taxonomy = ?)