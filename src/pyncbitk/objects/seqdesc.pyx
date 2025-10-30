# cython: language_level=3
"""Sequence description data model.

A `SeqDesc` is meant to describe a `BioSeq` in a biological and/or 
bibliographic context.

See Also:
    The `Seq-descr <https://ncbi.github.io/cxx-toolkit/pages/ch_datamod#ch_datamod.Seqdescr_Describing_>`_
    section of the NCBI C++ Toolkit documentation.

"""

from libcpp.algorithm cimport swap
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.list cimport list as cpplist

from cpython cimport Py_buffer
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.buffer cimport PyBUF_FORMAT
from cpython.object cimport Py_LT, Py_EQ, Py_GT, Py_LE, Py_NE, Py_GE

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seqdesc cimport CSeqdesc, E_Choice as CSeqdesc_choice
from ..toolkit.serial.serialbase cimport CSerialObject

from ..serial cimport Serial

import pickle


# --- SeqDesc ------------------------------------------------------------------

cdef class SeqDesc(Serial):
    """An abstract base storage of sequence description.
    """

    @staticmethod
    cdef SeqDesc _wrap(CRef[CSeqdesc] ref):
        cdef SeqDesc         obj
        cdef CSeqdesc_choice kind = ref.GetNonNullPointer().Which()

        if kind == CSeqdesc_choice.e_not_set:
            obj = SeqDesc.__new__(SeqDesc)
        elif kind == CSeqdesc_choice.e_Name:
            obj = NameDesc.__new__(NameDesc)
        elif kind == CSeqdesc_choice.e_Title:
            obj = TitleDesc.__new__(TitleDesc)
        elif kind == CSeqdesc_choice.e_Region:
            obj = RegionDesc.__new__(RegionDesc)
        else:
            obj = SeqDesc.__new__(SeqDesc)
            # return NotImplemented

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(self):
        self._ref.Reset(new CSeqdesc())
        self._ref.GetNonNullPointer().Select(CSeqdesc_choice.e_not_set)


cdef class NameDesc(SeqDesc):
    """A description storing the name of a sequence.
    """

    def __init__(self, str data):
        cdef bytes                    _data
        cdef const unsigned char[::1] _view
        cdef string                   s

        if isinstance(data, str):
            _data = data.encode()
            _view = _data
        else:
            _view = data

        super().__init__()

        with nogil:
            s = string(<const char*> &_view[0], _view.shape[0])
            self._ref.GetNonNullPointer().Select(CSeqdesc_choice.e_Name)
            self._ref.GetNonNullPointer().SetName(s)

    def __reduce__(self):
        return type(self), (str(self),)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({str(self)!r})"

    def __rich_repr__(self):
        yield str(self)

    def __str__(self):
        return self._ref.GetNonNullPointer().GetName().decode()


cdef class TitleDesc(SeqDesc):
    """A description storing the title of a sequence.
    """

    def __init__(self, str data):
        cdef bytes                    _data
        cdef const unsigned char[::1] _view
        cdef string                   s

        if isinstance(data, str):
            _data = data.encode()
            _view = _data
        else:
            _view = data

        super().__init__()

        with nogil:
            s = string(<const char*> &_view[0], _view.shape[0])
            self._ref.GetNonNullPointer().Select(CSeqdesc_choice.e_Title)
            self._ref.GetNonNullPointer().SetTitle(s)

    def __reduce__(self):
        return type(self), (str(self),)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({str(self)!r})"

    def __rich_repr__(self):
        yield str(self)

    def __str__(self):
        return self._ref.GetNonNullPointer().GetTitle().decode()


cdef class RegionDesc(SeqDesc):

    def __init__(self, str data):
        cdef bytes                    _data
        cdef const unsigned char[::1] _view
        cdef string                   s

        if isinstance(data, str):
            _data = data.encode()
            _view = _data
        else:
            _view = data

        super().__init__()

        with nogil:
            s = string(<const char*> &_view[0], _view.shape[0])
            self._ref.GetNonNullPointer().Select(CSeqdesc_choice.e_Region)
            self._ref.GetNonNullPointer().SetRegion(s)

    def __reduce__(self):
        return type(self), (str(self),)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({str(self)!r})"

    def __rich_repr__(self):
        yield str(self)

    def __str__(self):
        return self._ref.GetNonNullPointer().GetRegion().decode()

# --- SeqDescSet ---------------------------------------------------------------

cdef class SeqDescSet(Serial):
    """A set of sequence descriptions.
    """
    
    @staticmethod
    cdef SeqDescSet _wrap(CRef[CSeq_descr] ref):
        cdef SeqDescSet obj = SeqDescSet.__new__(SeqDescSet)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(self, items=()):
        """__init__(self, items=())\n--\n

        Create a new set of sequence descriptions.

        Arguments:
            items (iterable of `~pyncbitk.objects.seq.SeqDesc`): An iterable
                of sequence descriptions to add to the set.

        """
        cdef SeqDesc item
        cdef cpplist[CRef[CSeqdesc]]* data

        self._ref.Reset(new CSeq_descr())
        data = &self._ref.GetNonNullPointer().GetMut()

        for item in items:
            data.push_back(item._ref)

    def __len__(self):
        return self._ref.GetNonNullPointer().Get().size()

    def __iter__(self):
        cdef cpplist[CRef[CSeqdesc]]* data = &self._ref.GetNonNullPointer().GetMut()
        for item in data[0]:
            yield SeqDesc._wrap(item)

    def __rich_repr__(self):
        yield list(self)