# cython: language_level=3, linetrace=True, binding=True

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
        cdef ObjectId obj
        cdef CObject_id_choice kind = ref.GetPointer().Which()

        if kind == CObject_id_choice.e_Id:
            obj = IntId.__new__(IntId)
        elif kind == CObject_id_choice.e_Str:
            obj = StrId.__new__(StrId)
        else:
            raise NotImplementedError

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(self):
        raise TypeError("Can't instantiate abstract class ObjectId")

cdef class StrId(ObjectId):
    """An object identifier which is stored as a C++ string.
    """

    def __init__(self, str id):
        cdef bytes _id = id.encode()
        cdef CObject_id* obj = new CObject_id()
        obj.Select(CObject_id_choice.e_Str)
        obj.SetStr(_id)
        self._ref = CRef[CObject_id](obj)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    @property
    def id(self):
        """`str`: The value of the string identifier.
        """
        return self._ref.GetNonNullPointer().GetStr().decode()

cdef class IntId(ObjectId):
    """An object identifier which is stored as an integer.
    """

    def __init__(self, int id):
        cdef CObject_id* obj = new CObject_id()
        obj.Select(CObject_id_choice.e_Id)
        obj.SetId(id)
        self._ref = CRef[CObject_id](obj)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    @property
    def id(self):
        """`int`: The value of the integer identifier.
        """
        return self._ref.GetNonNullPointer().GetId()
