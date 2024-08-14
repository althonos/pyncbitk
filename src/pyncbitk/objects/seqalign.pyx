# cython: language_level=3, linetrace=True, binding=True

from libc.math cimport NAN
from libcpp cimport bool
from libcpp.string cimport string

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.general.object_id cimport CObject_id
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id
from ..toolkit.objects.seqloc.na_strand cimport ENa_strand
from ..toolkit.objects.seqalign.dense_seg cimport CDense_seg
from ..toolkit.objects.seqalign.seq_align cimport CSeq_align, EScoreType, EType as CSeq_align_type
from ..toolkit.objects.seqalign.seq_align cimport C_Segs, E_Choice as C_Segs_choice
from ..toolkit.objects.seqalign.seq_align_set cimport CSeq_align_set
from ..toolkit.objects.seqalign.score cimport CScore, C_Value as CScore_value, E_Choice as CScore_value_choice
from ..toolkit.serial.serialbase cimport CSerialObject

from ..serial cimport Serial
from .general cimport ObjectId
from .seqloc cimport SeqId

# --- Constants ----------------------------------------------------------------

cdef dict _NA_STRAND_STR = {
    ENa_strand.eNa_strand_unknown: "unknown",
    ENa_strand.eNa_strand_plus: "plus",
    ENa_strand.eNa_strand_minus: "minus",
    ENa_strand.eNa_strand_both: "both",
    ENa_strand.eNa_strand_both_rev: "bothrev",
    ENa_strand.eNa_strand_other: "other",
}

cdef dict _NA_STRAND_ENUM = {
    v:k for k,v in _NA_STRAND_STR.items()
}

# --- SeqAlign -----------------------------------------------------------------

cdef class SeqAlignScore:
    # TODO: inherit

    @staticmethod
    cdef SeqAlignScore _wrap(CRef[CScore] ref):
        cdef SeqAlignScore score = SeqAlignScore.__new__(SeqAlignScore)
        score._ref = ref
        return score

    def __iter__(self):
        yield self.id
        yield self.value

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(id={self.id!r}, value={self.value!r})"

    @property
    def id(self):
        if not self._ref.GetObject().IsSetId():
            return None
        id_ = &self._ref.GetObject().GetIdMut()
        cref = CRef[CObject_id](id_)
        return ObjectId._wrap(cref)

    @property
    def value(self):
        value = &self._ref.GetObject().GetValueMut()
        kind = value.Which()
        if kind == CScore_value_choice.e_Int:
            return value.GetInt()
        elif kind == CScore_value_choice.e_Real:
            return value.GetReal()
        raise TypeError(f"Unknown value type: {kind!r}")


cdef class AlignRow:
    
    @property
    def start(self):
        cdef CSeq_align* obj = &self._ref.GetObject()
        return obj.GetSeqStart(self._row)

    @property
    def stop(self):
        cdef CSeq_align* obj = &self._ref.GetObject()
        return obj.GetSeqStop(self._row)

    @property
    def id(self):
        cdef CSeq_align*    obj = &self._ref.GetObject()
        cdef const CSeq_id* id_ = &obj.GetSeq_id(self._row)
        return SeqId._wrap(CRef[CSeq_id](<CSeq_id*> id_))


cdef class AlignSegments(Serial):
    
    @staticmethod
    cdef AlignSegments _wrap(CRef[C_Segs] ref):
        cdef AlignSegments obj
        cdef C_Segs_choice kind = ref.GetObject().Which() 

        if kind == C_Segs.E_Choice.e_Denseg:
            obj = DenseSegments.__new__(DenseSegments)
        elif kind == C_Segs.E_Choice.e_not_set:
            obj = AlignSegments.__new__(AlignSegments)
        else:
            raise NotImplementedError

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()


cdef class DenseSegments(AlignSegments):
    
    @property
    def data(self):
        cdef CRef[CDense_seg] ref
        cdef CDense_seg* seg = &self._ref.GetObject().GetDensegMut()
        return DenseSegmentsData._wrap(CRef[CDense_seg](seg))

cdef class DenseSegmentsData(Serial):

    @staticmethod
    cdef DenseSegmentsData _wrap(CRef[CDense_seg] ref):
        cdef DenseSegmentsData obj = DenseSegmentsData.__new__(DenseSegmentsData)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    @property
    def num_segments(self):
        cdef CDense_seg* obj = &self._ref.GetObject()
        return obj.GetNumseg()

    @property
    def ids(self):
        """`list` of `SeqId`: The identifiers of the sequences in the segment.
        """
        cdef CDense_seg* obj = &self._ref.GetObject()
        cdef list        ids = []

        for ref in obj.GetIdsMut(): # FIXME: const iteration
            ids.append(SeqId._wrap(ref))
        return ids

    @property
    def lengths(self):
        """`list` of `int`: The lengths of each segment.
        """
        # FIXME: make zero-copy
        cdef CDense_seg* obj  = &self._ref.GetObject()
        cdef list        lens = []

        for length in obj.GetLensMut(): # FIXME: const iteration
            lens.append(length)
        return lens

    @property
    def starts(self):
        """`list` of `int`: The start offsets for the sequences in each segment.
        """
        # FIXME: make zero-copy
        cdef CDense_seg* obj    = &self._ref.GetObject()
        cdef list        starts = []

        for start in obj.GetStartsMut(): # FIXME: const iteration
            starts.append(start)
        return starts

    @property
    def strands(self):
        """`list` of `str`, or `None`: The strand for each sequence, if any.
        """
        cdef ENa_strand  strand
        cdef CDense_seg* obj     = &self._ref.GetObject()
        cdef list        strands = []

        for strand in obj.GetStrandsMut():
            strands.append(_NA_STRAND_STR[strand])
        return strands

cdef class SeqAlign(Serial):

    @staticmethod
    cdef SeqAlign _wrap(CRef[CSeq_align] ref):
        cdef SeqAlign        obj 
        cdef CSeq_align_type ty  = ref.GetObject().GetType()
        if ty == CSeq_align_type.eType_not_set:
            obj = SeqAlign.__new__(SeqAlign)
        elif ty == CSeq_align_type.eType_global:
            obj = GlobalSeqAlign.__new__(GlobalSeqAlign)
        elif ty == CSeq_align_type.eType_diags:
            obj = DiagonalSeqAlign.__new__(DiagonalSeqAlign)
        elif ty == CSeq_align_type.eType_partial:
            obj = PartialSeqAlign.__new__(PartialSeqAlign)
        elif ty == CSeq_align_type.eType_disc:
            obj = DiscontinuousSeqAlign.__new__(DiscontinuousSeqAlign)
        else:
            raise RuntimeError("Unsupported `SeqAlign` type")
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __len__(self):
        cdef CSeq_align* obj = &self._ref.GetObject()
        if not obj.IsSetDim():
            return 0
        return obj.GetDim()

    def __getitem__(self, ssize_t index):
        cdef CSeq_align* obj    = &self._ref.GetObject()
        cdef ssize_t     length = 0
        cdef ssize_t     index_ = index

        if obj.IsSetDim():
            length = obj.GetDim()

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)

        cdef AlignRow row = AlignRow.__new__(AlignRow)
        row._ref = self._ref
        row._row = index
        return row

    @property
    def evalue(self):
        cdef double evalue = NAN
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_EValue, evalue):
            return None
        return evalue

    @property
    def scores(self):
        cdef CRef[CScore]  ref
        cdef SeqAlignScore score
        cdef list          scores = []

        if not self._ref.GetObject().IsSetScore():
            return None

        for ref in self._ref.GetObject().GetScoreMut():
            scores.append(SeqAlignScore._wrap(ref))

        return scores

    @property
    def segments(self):
        cdef CSeq_align*  obj = &self._ref.GetObject()
        cdef CRef[C_Segs] ref = CRef[C_Segs](&obj.GetSegsMut())
        return AlignSegments._wrap(ref)




cdef class GlobalSeqAlign(SeqAlign):
    pass

cdef class DiagonalSeqAlign(SeqAlign):
    pass

cdef class PartialSeqAlign(SeqAlign):
    pass

cdef class DiscontinuousSeqAlign(SeqAlign):
    pass

cdef class SeqAlignSet(Serial):
    @staticmethod
    cdef SeqAlignSet _wrap(CRef[CSeq_align_set] ref):
        cdef SeqAlignSet obj = SeqAlignSet.__new__(SeqAlignSet)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __iter__(self):
        cdef CRef[CSeq_align] ref
        for ref in self._ref.GetObject().GetMut():
            yield SeqAlign._wrap(ref)

    def __len__(self):
        return self._ref.GetObject().Get().size()