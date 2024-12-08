from ....corelib.ncbiobj cimport CObject, CRef, CConstRef

cdef extern from "algo/blast/igblast/igblast.hpp" namespace "ncbi::blast" nogil:

    cppclass CIgBlastOptions(CObject):
        pass