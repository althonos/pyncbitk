# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_data cimport CSeq_data, E_Choice as CSeq_data_choice

from ..serial cimport Serial

# --- SeqData ------------------------------------------------------------------

cdef class SeqData(Serial):
    cdef CRef[CSeq_data] _ref

    @staticmethod
    cdef SeqData _wrap(CRef[CSeq_data] ref)
    cdef bool _validate(self) except False

    cpdef SeqData copy(self, bool pack=*)
    cpdef SeqData complement(self, bool pack=*)
    cpdef SeqData reverse_complement(self, bool pack=*)

cdef class SeqAaData(SeqData):
    cpdef str decode(self)

cdef class SeqNaData(SeqData):
    cpdef str decode(self)

cdef class IupacNaData(SeqNaData):
    pass

cdef class IupacAaData(SeqAaData):
    pass

cdef class Ncbi2NaData(SeqNaData):
    pass

cdef class Ncbi4NaData(SeqNaData):
    pass

cdef class Ncbi8NaData(SeqNaData):
    pass

cdef class NcbiPNaData(SeqNaData):
    pass

cdef class Ncbi8AaData(SeqAaData):
    pass

cdef class NcbiEAaData(SeqAaData):
    pass

cdef class NcbiPAaData(SeqAaData):
    pass

cdef class NcbiStdAa(SeqAaData):
    pass

cdef class GapData(SeqData):
    pass


