from libcpp cimport bool

from ....corelib.ncbiobj cimport CObject, CRef
from .blast_options_handle cimport CBlastOptionsHandle
from .blast_options cimport EAPILocality, CBlastOptions

cdef extern from "algo/blast/api/blast_nucl_options.hpp" namespace "ncbi::blast" nogil:
    
    cppclass CBlastNucleotideOptionsHandle(CBlastOptionsHandle):
        CBlastNucleotideOptionsHandle()
        CBlastNucleotideOptionsHandle(EAPILocality locality)
        CBlastNucleotideOptionsHandle(CRef[CBlastOptions] opt)

        void SetDefaults()
        void SetTraditionalBlastnDefaults()
        void SetTraditionalMegablastDefaults()

        bool GetDustFiltering()
        void SetDustFiltering(bool val)

        int GetMatchReward()
        void SetMatchReward(int r)
        int GetMismatchPenalty()
        void SetMismatchPenalty(int p)
        const char* GetMatrixName()
        void SetMatrixName(const char* matrix)
        int GetGapOpeningCost()
        void SetGapOpeningCost(int g)
        int GetGapExtensionCost()
        void SetGapExtensionCost(int e)