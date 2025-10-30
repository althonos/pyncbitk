# cython: language_level=3
"""Actual sequence data of biological sequences.

The NCBI C++ Toolkit provides a unified API for encoding the actual data of
protein and nucleotide sequences. It supports textual representation,
ordinal encoding (replacing *A*, *C*, *G*, *T* with *1*, *2*, *3*, *4*), as
well as compressed bitmap representations.

See Also:
    The `Data Model <https://ncbi.github.io/cxx-toolkit/pages/ch_datamod#ch_datamod.Seqdata_Encoding_the>`_
    chapter of the NCBI C++ Toolkit documentation.

"""

from libcpp.algorithm cimport swap
from libcpp.string cimport string
from libcpp.vector cimport vector

from cpython cimport Py_buffer
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.buffer cimport PyBUF_FORMAT
from cpython.object cimport Py_LT, Py_EQ, Py_GT, Py_LE, Py_NE, Py_GE

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_data cimport CSeq_data, E_Choice as CSeq_data_choice
from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.seq.seqport_util cimport CSeqportUtil
from ..toolkit.objects.seq.iupacna cimport CIUPACna
from ..toolkit.objects.seq.iupacaa cimport CIUPACaa
from ..toolkit.objects.seq.ncbi2na cimport CNCBI2na
from ..toolkit.objects.seq.ncbi4na cimport CNCBI4na
from ..toolkit.objects.seq.ncbi8na cimport CNCBI8na
from ..toolkit.objects.seq.ncbieaa cimport CNCBIeaa

from ..serial cimport Serial

import pickle

# --- SeqData ------------------------------------------------------------------

cdef class SeqData(Serial):
    """An abstract base storage of sequence data.
    """

    @staticmethod
    cdef SeqData _wrap(CRef[CSeq_data] ref):
        cdef SeqData          obj
        cdef CSeq_data_choice kind = ref.GetNonNullPointer().Which()

        if kind == CSeq_data_choice.e_not_set:
            obj = SeqData.__new__(SeqData)
        elif kind == CSeq_data_choice.e_Iupacna:
            obj = IupacNaData.__new__(IupacNaData)
        elif kind == CSeq_data_choice.e_Iupacaa:
            obj = IupacAaData.__new__(IupacAaData)
        elif kind == CSeq_data_choice.e_Ncbi2na:
            obj = Ncbi2NaData.__new__(Ncbi2NaData)
        elif kind == CSeq_data_choice.e_Ncbi4na:
            obj = Ncbi4NaData.__new__(Ncbi4NaData)
        elif kind == CSeq_data_choice.e_Ncbi8na:
            obj = Ncbi8NaData.__new__(Ncbi8NaData)
        elif kind == CSeq_data_choice.e_Ncbipna:
            obj = NcbiPNaData.__new__(NcbiPNaData)
        elif kind == CSeq_data_choice.e_Ncbi8aa:
            obj = Ncbi8AaData.__new__(Ncbi8AaData)
        elif kind == CSeq_data_choice.e_Ncbieaa:
            obj = NcbiEAaData.__new__(NcbiEAaData)
        elif kind == CSeq_data_choice.e_Ncbipaa:
            obj = NcbiPAaData.__new__(NcbiPAaData)
        elif kind == CSeq_data_choice.e_Ncbistdaa:
            obj = NcbiStdAa.__new__(NcbiStdAa)
        elif kind == CSeq_data_choice.e_Gap:
            obj = GapData.__new__(GapData)
        else:
            raise NotImplementedError

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    cdef bool _validate(self) except False:
        if not CSeqportUtil.FastValidate(self._ref.GetObject()):
            raise ValueError("Invalid elements in sequence data")
        return True

    @classmethod
    def _variant(cls):
        raise NotImplementedError("SeqData._variant")

    def __init__(self):
        self._ref.Reset(new CSeq_data())
        self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_not_set)

    def __copy__(self):
        return self.copy()

    cpdef SeqData complement(self, bool pack=False):  # FIXME: move to SeqNaData?
        cdef CSeq_data* data = new CSeq_data()
        with nogil:
            CSeqportUtil.Complement(self._ref.GetObject(), data)
            if pack:
                CSeqportUtil.Pack(data)
        return SeqData._wrap(CRef[CSeq_data](data))

    cpdef SeqData reverse_complement(self, bool pack=False):  # FIXME: move to SeqNaData?
        cdef CSeq_data* data = new CSeq_data()
        with nogil:
            CSeqportUtil.ReverseComplement(self._ref.GetObject(), data)
            if pack:
                CSeqportUtil.Pack(data)
        return SeqData._wrap(CRef[CSeq_data](data))

    cpdef SeqData copy(self, bool pack=False):
        """Create a copy of the sequence data.

        Arguments
            pack (`bool`): Whether to perform additional packing of the
                data stored in the copy, allowing a more compact
                representation but potentially changing the `SeqData`
                subtype.

        """
        cdef CSeq_data* data = new CSeq_data()
        with nogil:
            CSeqportUtil.GetCopy(self._ref.GetObject(), data)
            if pack:
                CSeqportUtil.Pack(data)
        return SeqData._wrap(CRef[CSeq_data](data))


cdef class SeqAaData(SeqData):
    """An abstract base storage of amino-acid sequence data.
    """

    @classmethod
    def encode(cls, object data):
        """Encode the textual sequence to a compressed representation.

        Arguments:
            data (`str`, `bytes`, or buffer-like object): The ASCII
                nucleotide sequence to be encoded. Python strings, and
                any other object supporting the buffer protocol is
                supported.

        """
        cdef CSeq_data_choice   variant = cls._variant()
        cdef SeqAaData          obj     = cls(bytearray())
        cdef NcbiEAaData        in_     = NcbiEAaData(data)

        with nogil:
            CSeqportUtil.Convert(
                in_._ref.GetObject(),
                obj._ref.GetNonNullPointer(),
                variant
            )

        return obj

    cpdef str decode(self):
        """Decode the contents of the sequence data.

        Returns:
            `str`: The decoded sequence data, as a Python string of
            NCBI-extended amino-acid symbols.

        """
        cdef CSeq_data*       out
        cdef CSeq_data*       data = self._ref.GetNonNullPointer()
        cdef CSeq_data_choice kind = data.Which()

        try:
            out = new CSeq_data()
            with nogil:
                CSeqportUtil.Convert(data[0], out, CSeq_data_choice.e_Ncbieaa)
            return out.GetNcbieaa().Get().decode()
        finally:
            del out


cdef class SeqNaData(SeqData):
    """An abstract base storage of nucleotide sequence data.
    """

    @classmethod
    def encode(cls, object data):
        """Encode the textual sequence to a compressed representation.

        Arguments:
            data (`str`, `bytes`, or buffer-like object): The ASCII
                nucleotide sequence to be encoded. Python strings, and
                any other object supporting the buffer protocol is
                supported.

        Raises:
            ValueError: When the sequence data contains invalid characters
                that do not belong to the nucleotide alphabet.

        """
        cdef CSeq_data_choice   variant = cls._variant()
        cdef SeqNaData          obj     = cls(bytearray())
        cdef IupacNaData        in_     = IupacNaData(data)

        with nogil:
            CSeqportUtil.Convert(
                in_._ref.GetObject(),
                obj._ref.GetNonNullPointer(),
                variant
            )

        return obj

    cpdef str decode(self):
        """Decode the contents of the sequence data.

        Returns:
            `str`: The decoded sequence data, as a Python string of
            IUPAC nucleotide symbols.

        """
        cdef CSeq_data*       out
        cdef CSeq_data*       data = self._ref.GetNonNullPointer()
        cdef CSeq_data_choice kind = data.Which()

        try:
            out = new CSeq_data()
            with nogil:
                CSeqportUtil.Convert(data[0], out, CSeq_data_choice.e_Iupacna)
            return out.GetIupacna().Get().decode()
        finally:
            del out


cdef class IupacNaData(SeqNaData):
    """Nucleotide sequence data stored as a IUPAC nucleotide string.

    Example:
        >>> seqdata = IupacNaData("ATTAGCCATGCATA")
        >>> seqdata.length
        14
        >>> seqdata.data
        b'ATTAGCCATGCATA'

    """

    @classmethod
    def _variant(cls):
        return CSeq_data_choice.e_Iupacna

    @classmethod
    def encode(cls, object data):
        return IupacNaData(data)

    def __init__(self, object data not None):
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
            self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_Iupacna)
            self._ref.GetNonNullPointer().SetIupacna(CIUPACna(s))

        self._validate()

    def __reduce_ex__(self, protocol):
        if protocol >= 5:
            return type(self), (pickle.PickleBuffer(self),)
        cdef bytes data = self._ref.GetNonNullPointer().GetIupacna().Get()
        return type(self), (data,)

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __rich_repr__(self):
        yield self.data

    def __richcmp__(self, other, int op):
        cdef IupacNaData   _other
        cdef const string* s1
        cdef const string* s2
        cdef int           status

        if not isinstance(other, IupacNaData):
            return NotImplemented

        _other = other
        s1 = &self._ref.GetNonNullPointer().GetIupacna().Get()
        s2 = &_other._ref.GetNonNullPointer().GetIupacna().Get()

        with nogil:
            status = s1[0].compare(s2[0])
        if op == Py_EQ:
            return status == 0
        elif op == Py_NE:
            return status != 0
        elif op == Py_GT:
            return status > 0
        elif op == Py_LE:
            return status <= 0
        elif op == Py_LT:
            return status < 0
        elif op == Py_GE:
            return status >= 0

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const string* data = &self._ref.GetNonNullPointer().GetIupacna().Get()

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL

        buffer.buf = <void*> data.data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = data.size() * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def length(self):
        """`int`: The length of the sequence data.
        """
        return self._ref.GetNonNullPointer().GetIupacna().Get().size()

    @property
    def data(self):
        """`bytes`: The sequence data as ASCII bytes.
        """
        cdef const string* data = &self._ref.GetNonNullPointer().GetIupacna().Get()
        return PyBytes_FromStringAndSize(data.data(), data.length())

    cpdef str decode(self):
        return self._ref.GetNonNullPointer().GetIupacna().Get().decode()


cdef class IupacAaData(SeqAaData):
    """Nucleotide sequence data stored in a IUPAC-UBI amino-acid string.

    The IUPAC-IUB Commission on Biochemical Nomenclature defined a code
    of one-letter abbreviations for the 20 standard amino-acids, as well
    as undeterminate and unknown symbols.

    References:
        * IUPAC-IUB Commission on Biochemical Nomenclature.
          "A One-Letter Notation for Amino Acid Sequences" 1–3. (1968).
          *Journal of Biological Chemistry*, 243(13), 3557–3559.
          :doi:`10.1016/S0021-9258(19)34176-6`.

    """

    @classmethod
    def _variant(cls):
        return CSeq_data_choice.e_Iupacaa

    @classmethod
    def encode(cls, object data):
        return IupacAaData(data)

    def __init__(self, object data):
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
            self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_Iupacaa)
            self._ref.GetNonNullPointer().SetIupacaa(CIUPACaa(s))

        self._validate()

    def __reduce_ex__(self, protocol):
        if protocol >= 5:
            return type(self), (pickle.PickleBuffer(self),), None
        return type(self), (memoryview(self).tobytes(),), None

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __rich_repr__(self):
        yield self.data

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const string* data = &self._ref.GetNonNullPointer().GetIupacna().Get()

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL

        buffer.buf = <void*> data.data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = data.size() * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def length(self):
        """`int`: The length of the sequence data.
        """
        return self._ref.GetNonNullPointer().GetIupacna().Get().size()

    @property
    def data(self):
        """`bytes`: The sequence data as ASCII bytes.
        """
        cdef const string* data = &self._ref.GetNonNullPointer().GetIupacna().Get()
        return PyBytes_FromStringAndSize(data.data(), data.length())

    cpdef str decode(self):
        return self._ref.GetNonNullPointer().GetIupacaa().Get().decode()


cdef class Ncbi2NaData(SeqNaData):
    """Nucleotide sequence data stored with 2-bit encoding.

    A nucleic acid containing no ambiguous bases can be encoded using a
    two-bit encoding per base, representing one of the four nucleobases:
    ``A``, ``C``, ``G`` or ``T``. This encoding is the most compact for
    unambiguous sequences.

    Example:
        >>> seqdata = Ncbi2NaData.encode("ATTAGCCATGCATA")
        >>> seqdata.data
        b'<\\x94\\xe4\\xc0'

    """

    @classmethod
    def _variant(cls):
        return CSeq_data_choice.e_Ncbi2na

    def __init__(self, object data not None):
        cdef CNCBI2na                 raw
        cdef const unsigned char[::1] view = data
        cdef size_t                   l    = view.shape[0]
        cdef vector[char]             vec  = vector[char]()

        super().__init__()

        with nogil:
            if l > 0:
                vec.insert(vec.end(), &view[0], &view[l-1])
            raw = CNCBI2na(vec)
            self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_Ncbi2na)
            swap[CNCBI2na](self._ref.GetNonNullPointer().GetNcbi2naMut(), raw)

        self._validate()

    def __reduce_ex__(self, protocol):
        if protocol >= 5:
            return type(self), (pickle.PickleBuffer(self),), None
        return type(self), (memoryview(self).tobytes(),), None

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __rich_repr__(self):
        yield self.data

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const vector[char]* data = &self._ref.GetObject().GetNcbi2na().Get()

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL

        buffer.buf = <void*> data.data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = data.size() * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def data(self):
        """`bytes`: The sequence data in 2-bit encoded format.
        """
        cdef const vector[char]* data = &self._ref.GetObject().GetNcbi2na().Get()
        return PyBytes_FromStringAndSize(data.data(), data.size())


cdef class Ncbi4NaData(SeqNaData):
    """Nucleotide sequence data stored with 4-bit encoding.
    """

    @classmethod
    def _variant(cls):
        return CSeq_data_choice.e_Ncbi4na

    def __init__(self, object data not None):
        cdef CNCBI4na        raw
        cdef const char[::1] view = data
        cdef size_t          l    = view.shape[0]
        cdef vector[char]    vec  = vector[char]()

        super().__init__()

        with nogil:
            if l > 0:
                vec.insert(vec.end(), &view[0], &view[l-1])
            raw = CNCBI4na(vec)
            self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_Ncbi4na)
            swap[CNCBI4na](self._ref.GetNonNullPointer().GetNcbi4naMut(), raw)

        self._validate()

    def __reduce_ex__(self, protocol):
        if protocol >= 5:
            return type(self), (pickle.PickleBuffer(self),), None
        return type(self), (memoryview(self).tobytes(),), None

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __rich_repr__(self):
        yield self.data

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const vector[char]* data = &self._ref.GetNonNullPointer().GetNcbi4na().Get()

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL

        buffer.buf = <void*> data.data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = data.size() * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def data(self):
        """`bytes`: The sequence data in 4-bit encoded format.
        """
        cdef const vector[char]* data = &self._ref.GetNonNullPointer().GetNcbi4na().Get()
        return PyBytes_FromStringAndSize(data.data(), data.size())


cdef class Ncbi8NaData(SeqNaData):

    @classmethod
    def _variant(cls):
        return CSeq_data_choice.e_Ncbi8na

    def __init__(self, object data not None):
        cdef CNCBI8na        raw
        cdef const char[::1] view = data
        cdef size_t          l    = view.shape[0]
        cdef vector[char]    vec  = vector[char]()

        super().__init__()

        with nogil:
            if l > 0:
                vec.insert(vec.end(), &view[0], &view[l-1])
            raw = CNCBI8na(vec)
            self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_Ncbi8na)
            swap[CNCBI8na](self._ref.GetNonNullPointer().GetNcbi8naMut(), raw)

        self._validate()

    def __reduce_ex__(self, protocol):
        if protocol >= 5:
            return type(self), (pickle.PickleBuffer(self),)
        return type(self), (memoryview(self).tobytes(),)

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __rich_repr__(self):
        yield self.data

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const vector[char]* data = &self._ref.GetNonNullPointer().GetNcbi8na().Get()

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL

        buffer.buf = <void*> data.data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = data.size() * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def data(self):
        """`bytes`: The sequence data in 8-bit encoded format.
        """
        cdef const vector[char]* data = &self._ref.GetNonNullPointer().GetNcbi8na().Get()
        return PyBytes_FromStringAndSize(data.data(), data.size())


cdef class NcbiPNaData(SeqNaData):
    """Nucleotide sequence data storing probabilities for each position.
    """


cdef class Ncbi8AaData(SeqAaData):
    """Amino-acid sequence data with support for modified residues.
    """


cdef class NcbiEAaData(SeqAaData):
    """Amino-acid sequence data storing an NCBI-extended string.

    This representation adds symbols for the non-standard selenocysteine
    amino-acid (`U`) as well as support for termination or gap characters.

    """

    @classmethod
    def _variant(cls):
        return CSeq_data_choice.e_Ncbieaa

    @classmethod
    def encode(cls, object data):
        return NcbiEAaData(data)

    def __init__(self, object data not None):
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
            self._ref.GetNonNullPointer().Select(CSeq_data_choice.e_Ncbieaa)
            self._ref.GetNonNullPointer().SetNcbieaa(CNCBIeaa(s))

        self._validate()

    def __reduce_ex__(self, protocol):
        if protocol >= 5:
            return type(self), (pickle.PickleBuffer(self),), None
        return type(self), (memoryview(self).tobytes(),), None

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __rich_repr__(self):
        yield self.data

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const string* data = &self._ref.GetNonNullPointer().GetNcbieaa().Get()

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL

        buffer.buf = <void*> data.data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = data.size() * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = NULL
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def length(self):
        """`int`: The length of the sequence data.
        """
        return self._ref.GetNonNullPointer().GetNcbieaa().Get().size()

    @property
    def data(self):
        """`bytes`: The sequence data as ASCII bytes.
        """
        cdef const string* data = &self._ref.GetNonNullPointer().GetNcbieaa().Get()
        return PyBytes_FromStringAndSize(data.data(), data.length())

    cpdef str decode(self):
        return self._ref.GetNonNullPointer().GetNcbieaa().Get().decode()


cdef class NcbiPAaData(SeqAaData):
    """Amino-acid sequence data storing probabilities for each position.
    """


cdef class NcbiStdAa(SeqAaData):
    """Amino-acid sequence data stored as ordinal encoding.

    This encoding represents the NCBI-extended amino-acids as consecutive
    integer values, starting with *0* for the gap character.

    """


cdef class GapData(SeqData):
    """A virtual sequence data storage representing a gap.
    """

