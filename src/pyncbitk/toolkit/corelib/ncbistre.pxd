from iostream cimport istream, ostream

cdef extern from "corelib/ncbistr.hpp" namespace "ncbi" nogil:

    ctypedef istream CNcbiIstream
    ctypedef ostream CNcbiOstream

