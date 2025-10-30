# cython: language_level=3
"""General objects for the NCBI C++ object model.
"""
cimport cpython.object

from libcpp.string cimport string

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.ncbimisc cimport TTaxId
from ..toolkit.objects.general.object_id cimport CObject_id, E_Choice as CObject_id_choice
from ..toolkit.objects.general.dbtag cimport CDbtag
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

    def __rich_repr__(self):
        yield self.value

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

# --- DBTag --------------------------------------------------------------------

cdef class DBTag(Serial):
    """A database cross-reference.
    """

    @staticmethod
    cdef DBTag _wrap(CRef[CDbtag] ref):
        cdef DBTag obj = DBTag.__new__(DBTag)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(self, str db not None, object tag not None):
        """Create a new `DBTag` object.

        Arguments:
            db (`str`): The name of the database or system.
            tag (`~pyncbitk.objects.general.ObjectId`, `str` or `int`): The 
                identifier of the resource in the database.

            
        """
        cdef ObjectId oid
        cdef CDbtag*  dbtag
        cdef bytes    b

        if isinstance(tag, (int, str)):
            oid = ObjectId(tag)
        elif isinstance(tag, ObjectId):
            oid = tag
        else:
            ty = type(tag).__name__
            raise TypeError(f"expected str, int or ObjectId, found {ty}")

        b = db.encode()

        self._ref.Reset(new CDbtag())
        dbtag = self._ref.GetNonNullPointer()
        dbtag.SetDb(b)
        dbtag.SetTag(oid._ref.GetObject())

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.db!r}, {self.tag!r})"

    def __rich_repr__(self):
        yield self.db
        yield self.tag

    def __richcmp__(self, object other, int op):
        cdef DBTag other_
        cdef int res
        if isinstance(other, DBTag):
            other_ = other
            res = self._ref.GetObject().Compare(other_._ref.GetObject())
            if op == cpython.object.Py_LT:
                return res < 0
            elif op == cpython.object.Py_LE:
                return res <= 0
            elif op == cpython.object.Py_EQ:
                return res == 0
            elif op == cpython.object.Py_GT:
                return res > 0
            elif op == cpython.object.Py_GE:
                return res >= 0
            elif op == cpython.object.Py_NE:
                return res != 0
            else:
                raise ValueError(op)
        return NotImplemented

    @property
    def db(self):
        """`str`: The name of database or system.
        """
        cdef string db = self._ref.GetNonNullPointer().GetDb()
        return db.decode()

    @property
    def tag(self):
        """`~pyncbitk.objects.general.ObjectId`: The database tag.
        """
        return ObjectId._wrap(CRef[CObject_id](&self._ref.GetNonNullPointer().GetTagMut()))


    cpdef str url(self, object taxonomy = None):
        """Get a URL to the resource, if available.

        Arguments:
            taxonomy (`int` or `None`): An optional taxonomy ID to use
                for URL generation.

        Returns:
            `str` or `None`: A URL to the resource pointed by this
            database tag, if available.

        """
        cdef string url
        if taxonomy is None:
            url = self._ref.GetObject().GetUrl()
        elif isinstance(taxonomy, str):
            raise NotImplementedError("DBTag.url")
        elif isinstance(taxonomy, int):
            url = self._ref.GetObject().GetUrl(<TTaxId> taxonomy)
        else:
            ty = type(taxonomy).__name__
            raise TypeError(f"expected int, str or None, found {ty!r}")
        return None if url.empty() else url.decode()