
from ..corelib.ncbiobj cimport CObject
from .seq_vector_ci cimport CSeqVectorTypes


cdef extern from "objmgr/seq_vector.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeqVector(CObject, CSeqVectorTypes):
        CSeqVector()

    