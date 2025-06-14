from .blast_options_handle cimport CBlastOptionsHandle
from .blast_options cimport EAPILocality

cdef extern from "algo/blast/api/blast_prot_options.hpp" namespace "ncbi::blast" nogil:
    
    cppclass CBlastProteinOptionsHandle(CBlastOptionsHandle):
        CBlastProteinOptionsHandle()
        CBlastProteinOptionsHandle(EAPILocality locality)
        # CBlastNucleotideOptionsHandle(CRef[CBlastOptions] opt)

        double GetWordThreshold() noexcept
        void SetWordThreshold(double wt)
        int GetWordSize() noexcept
        void SetWordSize(int ws)

        void SetDefaults()
        void SetTraditionalBlastnDefaults()
        void SetTraditionalMegablastDefaults()