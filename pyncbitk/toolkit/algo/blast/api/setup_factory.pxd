from libcpp cimport bool

from ....corelib.ncbiobj cimport CObject, CRef

cdef extern from "algo/blast/api/setup_factory.hpp" namespace "ncbi::blast" nogil:


    cppclass CStructWrapper[TData](CObject):
        TData* GetPointer()

    cppclass CThreadable:
        CThreadable()

        void SetNumberOfThreads(size_t nthreads)
        size_t GetNumberOfThreads() const
        bool IsMultiThreaded() const