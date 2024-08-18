from libcpp cimport bool

from .blast_prot_options cimport CBlastProteinOptionsHandle
from .blast_options cimport EAPILocality

cdef extern from "algo/blast/api/blast_advprot_options.hpp" namespace "ncbi::blast" nogil:
    
    cppclass CBlastAdvancedProteinOptionsHandle(CBlastProteinOptionsHandle):
        CBlastAdvancedProteinOptionsHandle()
        CBlastAdvancedProteinOptionsHandle(EAPILocality locality)
        # CBlastNucleotideOptionsHandle(CRef[CBlastOptions] opt)

        # ECompoAdjustModes GetCompositionBasedStats() const
        # void SetCompositionBasedStats(ECompoAdjustModes mode)
        
        bool GetSmithWatermanMode() const
        void SetSmithWatermanMode()
        void SetSmithWatermanMode(bool m)