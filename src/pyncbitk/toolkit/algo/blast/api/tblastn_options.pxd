from libcpp cimport bool

from .blast_advprot_options cimport CBlastAdvancedProteinOptionsHandle
from .blast_options cimport EAPILocality

cdef extern from "algo/blast/api/tblastn_options.hpp" namespace "ncbi::blast" nogil:
    
    cppclass CTBlastnOptionsHandle(CBlastAdvancedProteinOptionsHandle):
        CTBlastnOptionsHandle()
        CTBlastnOptionsHandle(EAPILocality locality)
        # CTBlastnOptionsHandle(CRef[CBlastOptions] opt)

        bool GetOutOfFrameMode()
        void SetOutOfFrameMode() noexcept
        void SetOutOfFrameMode(bool m) noexcept

        int GetFrameShiftPenalty() const 
        void SetFrameShiftPenalty(int p) except +

        int GetLongestIntronLength() const
        void SetLongestIntronLength(int l) except +

        int GetDbGeneticCode() const
        void SetDbGeneticCode(int gc) except +