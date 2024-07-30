from .blast_options_handle cimport CBlastOptionsHandle
from .blast_options cimport EAPILocality

cdef extern from "algo/blast/api/blast_nucl_options.hpp" namespace "ncbi::blast" nogil:
    
    cppclass CBlastNucleotideOptionsHandle(CBlastOptionsHandle):
        CBlastNucleotideOptionsHandle()
        CBlastNucleotideOptionsHandle(EAPILocality locality)
        # CBlastNucleotideOptionsHandle(CRef[CBlastOptions] opt)

        void SetDefaults()
        void SetTraditionalBlastnDefaults()
        void SetTraditionalMegablastDefaults()