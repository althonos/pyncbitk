from ..corelib.ncbiobj cimport CObject

cdef extern from "objmgr/data_loader.hpp" namespace "ncbi::objects" nogil:

    cppclass CDataLoader(CObject):
        CDataLoader()