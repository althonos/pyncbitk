from libcpp cimport bool
from libcpp.string cimport string

cdef extern from * namespace "std" nogil:
    ctypedef ssize_t streamsize
    ctypedef ssize_t streampos
    cdef cppclass ios_base:
        cppclass openmode:
            openmode operator|(openmode)

cdef extern from * namespace "std" nogil:
    cdef cppclass istream:
        istream(streambuf* sb)
        istream& seekg(streampos pos)

cdef extern from * namespace "std" nogil:
    cdef cppclass streambuf:
        char* eback()
        char* egptr()
        streambuf* pubsetbuf(char* s, streamsize n);
        streambuf* setbuf(char* s, streamsize n);
        void setg(char* gbeg, char* gnext, char* gend);

cdef extern from * namespace "std" nogil:
    cdef cppclass filebuf(streambuf):
        filebuf()
        filebuf* open(const char* filename, ios_base.openmode mode)
        bool is_open()
        filebuf* close()

cdef extern from * namespace "std" nogil:
    cdef cppclass stringbuf(streambuf):
        stringbuf()
        stringbuf(const string& str)
        string str() const
        void str (const string& str)

cdef extern from * namespace "std" nogil:
    cdef cppclass ostream:
        ostream(streambuf* sb)

cdef extern from "corelib/ncbistr.hpp" namespace "ncbi" nogil:

    ctypedef istream CNcbiIstream
    ctypedef ostream CNcbiOstream

