# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.general.object_id cimport CObject_id, E_Choice as CObject_id_choice
from ..toolkit.serial.serialbase cimport CSerialObject

from ..serial cimport Serial

# --- ObjectId -----------------------------------------------------------------

ctypedef fused ObjectIdValue:
    int
    str
    bytes

cdef class ObjectId(Serial):
    """A basic identifier for any NCBI Toolkit object.
    """

    @staticmethod
    cdef ObjectId _wrap(CRef[CObject_id] ref):
        cdef ObjectId obj = ObjectId.__new__(ObjectId)
        obj._ref = ref
        # cdef CObject_id_choice kind = ref.GetPointer().Which()
        # if kind == CObject_id_choice.e_Id:
        #     obj = IntId.__new__(IntId)
        # elif kind == CObject_id_choice.e_Str:
        #     obj = StrId.__new__(StrId)
        # else:
        #     raise NotImplementedError
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(self, ObjectId value):
        cdef CObject_id* obj = new CObject_id()
        if ObjectId is int:
            obj.Select(CObject_id_choice.e_Id)
            obj.SetId(value)
        elif ObjectId is str:
            obj.Select(CObject_id_choice.e_Str)
            obj.SetStr(value.encode())
        elif ObjectId is bytes:
            obj.Select(CObject_id_choice.e_Str)
            obj.SetStr(value)
        self._ref.Reset(obj)

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.value!r})"

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

# cdef class StrId(ObjectId):
#     """An object identifier which is stored as a C++ string.
#     """

#     def __init__(self, str id):
#         cdef bytes _id = id.encode()
#         cdef CObject_id* obj = new CObject_id()
#         obj.Select(CObject_id_choice.e_Str)
#         obj.SetStr(_id)
#         self._ref = CRef[CObject_id](obj)

#     def __repr__(self):
#         cdef str ty = type(self).__name__
#         return f"{ty}({self.id!r})"

#     @property
#     def id(self):
#         """`str`: The value of the string identifier.
#         """
#         return self._ref.GetNonNullPointer().GetStr().decode()

# cdef class IntId(ObjectId):
#     """An object identifier which is stored as an integer.
#     """

#     def __init__(self, int id):
#         cdef CObject_id* obj = new CObject_id()
#         obj.Select(CObject_id_choice.e_Id)
#         obj.SetId(id)
#         self._ref = CRef[CObject_id](obj)

#     def __repr__(self):
#         cdef str ty = type(self).__name__
#         return f"{ty}({self.id!r})"

#     @property
#     def id(self):
#         """`int`: The value of the integer identifier.
#         """
#         return self._ref.GetNonNullPointer().GetId()
