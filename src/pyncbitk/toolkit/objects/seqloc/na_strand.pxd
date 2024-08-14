from libcpp cimport bool

cdef extern from "objects/seqloc/Na_strand_.hpp" namespace "ncbi::objects" nogil:

    enum ENa_strand:
        eNa_strand_unknown 
        eNa_strand_plus    
        eNa_strand_minus   
        eNa_strand_both     
        eNa_strand_both_rev 
        eNa_strand_other   

cdef extern from "objects/seqloc/Na_strand.hpp" namespace "ncbi::objects" nogil:
    
    enum EIsSetStrand:
        eIsSetStrand_Any
        eIsSetStrand_All

    enum ESeqLocExtremes:
        eExtreme_Biological
        eExtreme_Positional

    bool IsForward(ENa_strand s) noexcept
    bool IsReverse(ENa_strand s) noexcept
    bool SameOrientation(ENa_strand a, ENa_strand b) noexcept
    ENa_strand Reverse(ENa_strand s) noexcept
   