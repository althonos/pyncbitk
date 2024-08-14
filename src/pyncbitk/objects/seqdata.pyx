# cython: language_level=3, linetrace=True, binding=True

from libcpp.string cimport string
from libcpp.vector cimport vector

from cpython cimport Py_buffer
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.buffer cimport PyBUF_FORMAT

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_data cimport CSeq_data, E_Choice as CSeq_data_choice
from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.seq.seqport_util cimport CSeqportUtil
from ..toolkit.objects.seq.iupacna cimport CIUPACna

from ..serial cimport Serial

# --- SeqData ------------------------------------------------------------------

cdef class SeqData(Serial):

    @staticmethod
    cdef SeqData _wrap(CRef[CSeq_data] ref):
        cdef SeqData          obj
        cdef CSeq_data_choice kind = ref.GetObject().Which()

        if kind == CSeq_data_choice.e_Iupacna:
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

    def __init__(self):
        self._ref.Reset(new CSeq_data())

    cpdef SeqData complement(self, bool pack=False):
        cdef CSeq_data* data = new CSeq_data()
        with nogil:
            CSeqportUtil.Complement(self._ref.GetObject(), data)
            if pack:
                CSeqportUtil.Pack(data)
        return SeqData._wrap(CRef[CSeq_data](data))

    cpdef SeqData reverse_complement(self, bool pack=False):
        cdef CSeq_data* data = new CSeq_data()
        with nogil:
            CSeqportUtil.ReverseComplement(self._ref.GetObject(), data)
            if pack:
                CSeqportUtil.Pack(data)
        return SeqData._wrap(CRef[CSeq_data](data))

    cpdef SeqData copy(self, bool pack=False):
        cdef CSeq_data* data = new CSeq_data()
        with nogil:
            CSeqportUtil.GetCopy(self._ref.GetObject(), data)
            if pack:
                CSeqportUtil.Pack(data)
        return SeqData._wrap(CRef[CSeq_data](data))


cdef class SeqAaData(SeqData):
    
    def decode(self):
        cdef CSeq_data*       out  
        cdef CSeq_data*       data = &self._ref.GetObject()
        cdef CSeq_data_choice kind = data.Which()

        try:
            out = new CSeq_data()
            with nogil:
                CSeqportUtil.Convert(data[0], out, CSeq_data_choice.e_Ncbieaa)
            return out.GetNcbieaa().Get().decode()
        finally:
            del out

cdef class SeqNaData(SeqData):
    
    def decode(self):
        cdef CSeq_data*       out  
        cdef CSeq_data*       data = &self._ref.GetObject()
        cdef CSeq_data_choice kind = data.Which()

        try:
            out = new CSeq_data()
            with nogil:
                CSeqportUtil.Convert(data[0], out, CSeq_data_choice.e_Iupacna)
            return out.GetIupacna().Get().decode()
        finally:
            del out

cdef class IupacNaData(SeqNaData):

    def __init__(self, object data):
        cdef bytes      _data

        if isinstance(data, str):
            _data = data.encode()
        else:
            _data = data

        super().__init__()
        self._ref.GetObject().Select(CSeq_data_choice.e_Iupacna)
        self._ref.GetObject().SetIupacna(CIUPACna(<string> _data))

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.data!r})"

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const string* data = &self._ref.GetObject().GetIupacna().Get()

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
        return self._ref.GetObject().GetIupacna().Get().size()

    @property
    def data(self):
        return self.decode()

    def decode(self):
        return self._ref.GetObject().GetIupacna().Get().decode()


cdef class IupacAaData(SeqAaData):
    
    def decode(self):
        return self._ref.GetObject().GetIupacaa().Get().decode()

cdef class Ncbi2NaData(SeqNaData):
    pass

cdef class Ncbi4NaData(SeqNaData):

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        cdef const vector[char]* data = &self._ref.GetObject().GetNcbi4na().Get()

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
        cdef const vector[char]* data = &self._ref.GetObject().GetNcbi4na().Get()
        return PyBytes_FromStringAndSize(data.data(), data.size())


cdef class Ncbi8NaData(SeqNaData):
    pass

cdef class NcbiPNaData(SeqNaData):
    pass

cdef class Ncbi8AaData(SeqAaData):
    pass

cdef class NcbiEAaData(SeqAaData):
    def decode(self):
        return self._ref.GetObject().GetNcbieaa().Get().decode()

cdef class NcbiPAaData(SeqAaData):
    pass

cdef class NcbiStdAa(SeqAaData):
    pass

cdef class GapData(SeqData):
    pass

