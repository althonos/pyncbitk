

cdef extern from "objmgr/seq_vector_ci.hpp" namespace "ncbi::objects::CSeqVectorTypes" nogil:

    ctypedef unsigned char       TResidue
    # ctypedef CSeq_data::E_Choice TCoding
    # ctypedef TResidue            value_type
    # ctypedef TSeqPos             size_type
    # ctypedef TSignedSeqPos       difference_type
    # ctypedef std::random_access_iterator_tag iterator_category
    ctypedef const TResidue*     pointer
    ctypedef const TResidue&     reference

cdef extern from "objmgr/seq_vector_ci.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeqVectorTypes:
        pass

    