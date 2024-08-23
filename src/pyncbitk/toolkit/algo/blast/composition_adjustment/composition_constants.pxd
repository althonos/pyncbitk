

cdef extern from "algo/blast/composition_adjustment/composition_constants.h" nogil:
    
    enum ECompoAdjustModes:
        eNoCompositionBasedStats
        eCompositionBasedStats
        eCompositionMatrixAdjust
        eCompoForceFullMatrixAdjust
        eNumCompoAdjustModes

    enum EMatrixAdjustRule:
        eDontAdjustMatrix             
        eCompoScaleOldMatrix           
        eUnconstrainedRelEntropy      
        eRelEntropyOldMatrixNewContext
        eRelEntropyOldMatrixOldContext
        eUserSpecifiedRelEntropy