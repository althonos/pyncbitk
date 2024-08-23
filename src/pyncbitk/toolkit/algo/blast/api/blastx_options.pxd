from libcpp cimport bool

from ....objects.seqloc.na_strand cimport ENa_strand
from ..composition_adjustment.composition_constants cimport ECompoAdjustModes
from .blast_prot_options cimport CBlastProteinOptionsHandle
from .blast_options cimport EAPILocality

cdef extern from "algo/blast/api/blastx_options.hpp" namespace "ncbi::blast" nogil:
    
    cppclass CBlastxOptionsHandle(CBlastProteinOptionsHandle):
        CBlastxOptionsHandle()
        CBlastxOptionsHandle(EAPILocality locality)
        # CTBlastnOptionsHandle(CRef[CBlastOptions] opt)

        ENa_strand GetStrandOption()
        void SetStrandOption(ENa_strand strand)

        int GetQueryGeneticCode() const
        void SetQueryGeneticCode(int gc)

        bool GetOutOfFrameMode()
        void SetOutOfFrameMode(bool m)        

        int GetFrameShiftPenalty()
        void SetFrameShiftPenalty(int p)

        int GetLongestIntronLength()
        void SetLongestIntronLength(int l)

        ECompoAdjustModes GetCompositionBasedStats()
        void SetCompositionBasedStats(ECompoAdjustModes mode)

        bool GetSmithWatermanMode()
        void SetSmithWatermanMode(bool m = false)

