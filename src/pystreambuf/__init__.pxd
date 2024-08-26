from iostream cimport streambuf

cdef extern from "fileobj/pywritebuf.h" nogil:

    cdef cppclass pywritebuf(streambuf):
        pywritebuf(object)


cdef extern from "fileobj/pyreadbuf.h" nogil:

    cdef cppclass pyreadbuf(streambuf):
        pyreadbuf(object)


cdef extern from "fileobj/pyreadintobuf.h" nogil:

    cdef cppclass pyreadintobuf(streambuf):
        pyreadintobuf(object)
