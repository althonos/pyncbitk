from libcpp cimport bool
from libcpp.vector cimport vector

from ...serial.serialbase cimport CSerialObject
from ...corelib.ncbiobj cimport CRef
from ...corelib.ncbimisc cimport TSeqPos, TSignedSeqPos
from ..seqloc.na_strand cimport ENa_strand
from ..seqloc.seq_id cimport CSeq_id
from .score cimport CScore

cdef extern from "objects/seqalign/Dense_seg_.hpp" namespace "ncbi::objects::CDense_seg_Base" nogil:

    ctypedef int TDim
    ctypedef int TNumseg
    ctypedef vector[CRef[CSeq_id]] TIds
    ctypedef vector[TSignedSeqPos] TStarts
    ctypedef vector[TSeqPos] TLens
    ctypedef vector[ENa_strand] TStrands
    ctypedef vector[CRef[CScore]] TScores



cdef extern from "objects/seqalign/Dense_seg_.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CDense_seg_Base(CSerialObject):
        CDense_seg_Base()

        bool IsSetDim() const
        bool CanGetDim() const
        void ResetDim()
        void SetDefaultDim()
        TDim GetDim() const
        void SetDim(TDim value)
        TDim& GetDimMut "SetDim"()

        bool IsSetNumseg() const
        bool CanGetNumseg() const
        void ResetNumseg()
        TNumseg GetNumseg() const
        void SetNumseg(TNumseg value)
        TNumseg& GetNumseg "SetNumseg"()

        bool IsSetIds() const
        bool CanGetIds() const
        void ResetIds()
        const TIds& GetIds() const
        TIds& GetIdsMut "SetIds"()

        bool IsSetStarts() const
        bool CanGetStarts() const
        void ResetStarts()
        const TStarts& GetStarts() const
        TStarts& GetStartsMut "SetStarts"()

        bool IsSetLens() const
        bool CanGetLens() const
        void ResetLens()
        const TLens& GetLens() const
        TLens& GetLensMut "SetLens"()

        bool IsSetStrands() const
        bool CanGetStrands() const
        void ResetStrands()
        const TStrands& GetStrands() const
        TStrands& GetStrandsMut "SetStrands"()

        bool IsSetScores() const
        bool CanGetScores() const
        void ResetScores()
        const TScores& GetScores() const
        TScores& GetScoresMut "SetScores"()

        void Reset()


cdef extern from "objects/seqalign/Dense_seg.hpp" namespace "ncbi::objects" nogil:

    cppclass CDense_seg(CDense_seg_Base):
        CDense_seg()