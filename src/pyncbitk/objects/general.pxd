from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.general.object_id cimport CObject_id

from ..serial cimport Serial


# --- ObjectId -----------------------------------------------------------------

cdef class ObjectId(Serial):
    cdef CRef[CObject_id] _ref

    @staticmethod
    cdef ObjectId _wrap(CRef[CObject_id] ref)
     
cdef class StrId(ObjectId):
    pass

cdef class IntId(ObjectId):
    pass