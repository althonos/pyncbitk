# cython: language_level=3
"""General objects for the NCBI C++ object model.
"""

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.general.object_id cimport CObject_id, E_Choice as CObject_id_choice
from ..toolkit.serial.serialbase cimport CSerialObject

from ..serial cimport Serial

# --- ObjectId -----------------------------------------------------------------

cdef class ObjectId(Serial):
    """A basic identifier for any NCBI Toolkit object.
    """

    @staticmethod
    cdef ObjectId _wrap(CRef[CObject_id] ref):
        cdef ObjectId obj = ObjectId.__new__(ObjectId)
        obj._ref = ref
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(self, object value):
        cdef bytes       _b
        cdef CObject_id* obj

        self._ref.Reset(new CObject_id())
        obj = self._ref.GetNonNullPointer()

        if isinstance(value, int):
            obj.Select(CObject_id_choice.e_Id)
            obj.SetId(value)
        elif isinstance(value, str):
            _b = value.encode()
            obj.Select(CObject_id_choice.e_Str)
            obj.SetStr(_b)
        else:
            _b = value
            obj.Select(CObject_id_choice.e_Str)
            obj.SetStr(_b)

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.value!r})"

    def __reduce__(self):
        return type(self), (self.value,)

    def __str__(self):
        return str(self.value)

    @property
    def value(self):
        """`str` or `int`: The actual value of the object identifier.
        """
        cdef CObject_id*       obj  = self._ref.GetNonNullPointer()
        cdef CObject_id_choice kind = obj.Which()
        if kind == CObject_id_choice.e_Id:
            return obj.GetId()
        else:
            return obj.GetStr().decode()
