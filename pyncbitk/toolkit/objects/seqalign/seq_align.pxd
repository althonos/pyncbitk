from libcpp cimport bool
from libcpp.list cimport list
from libcpp.vector cimport vector

from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from ..general.object_id cimport CObject_id
from ..seqloc.seq_loc cimport CSeq_loc
from .score cimport CScore


# cdef extern from "objects/seqalign/Seq_align_.hpp" namespace "ncbi::objects::CSeq_align_Base::C_Segs" nogil:

#     ctypedef CSerialObject Tparent



cdef extern from "objects/seqalign/Seq_align_.hpp" namespace "ncbi::objects::CSeq_align_Base" nogil:

    ctypedef CSerialObject Tparent
    
    enum EType:
        eType_not_set
        eType_global
        eType_diags
        eType_partial
        eType_disc
        eType_other

    cppclass C_Segs(CSerialObject):
        pass

    ctypedef EType TType
    ctypedef int TDim
    ctypedef vector[CRef[CScore]] TScore
    ctypedef C_Segs TSegs
    ctypedef list[CRef[CSeq_loc]] TBounds
    ctypedef list[CRef[CObject_id]] TId
    # ctypedef list[CRef[CUser_object]] TExt

    enum E_memberIndex:
        e__allMandatory
        e_type
        e_dim
        e_score
        e_segs
        e_bounds
        e_id
        e_ext

    

cdef extern from "objects/seqalign/Seq_align_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_align_Base(CSerialObject):
        CSeq_align_Base()

        bool IsSetType() const
        bool CanGetType() const
        void ResetType()
        TType GetType() const
        void SetType(TType value)
        TType& GetTypeRw "SetType" ()

        bool IsSetDim() const
        bool CanGetDim() const
        void ResetDim()
        TDim GetDim() const
        void SetDim(TDim value)
        TDim& GetDimRw "SetDim" ()

        bool IsSetScore() const
        bool CanGetScore() const
        void ResetScore()
        const TScore& GetScore() const
        TScore& GetScoreRw "SetScore" ()


cdef extern from "objects/seqalign/Seq_align.hpp" namespace "ncbi::objects::CSeq_align" nogil:

    ctypedef CSeq_align_Base Tparent
    ctypedef int TDim


cdef extern from "objects/seqalign/Seq_align.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_align(CSeq_align_Base):
        CSeq_align()
