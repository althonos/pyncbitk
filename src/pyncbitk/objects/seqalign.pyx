# cython: language_level=3, linetrace=True, binding=True

from libc.math cimport NAN
from libcpp cimport bool
from libcpp.string cimport string

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.ncbimisc cimport TSeqPos
from ..toolkit.objects.general.object_id cimport CObject_id
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id
from ..toolkit.objects.seqloc.na_strand cimport ENa_strand
from ..toolkit.objects.seqalign.dense_seg cimport CDense_seg
from ..toolkit.objects.seqalign.seq_align cimport CSeq_align, EScoreType, EType as CSeq_align_type
from ..toolkit.objects.seqalign.seq_align cimport C_Segs, E_Choice as C_Segs_choice
from ..toolkit.objects.seqalign.seq_align_set cimport CSeq_align_set
from ..toolkit.objects.seqalign.score cimport CScore, C_Value as CScore_value, E_Choice as CScore_value_choice
from ..toolkit.objtools.align_format.align_format_util cimport CAlignFormatUtil
from ..toolkit.serial.serialbase cimport CSerialObject

from ..serial cimport Serial
from .general cimport ObjectId
from .seqid cimport SeqId

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


# --- Accessory classes --------------------------------------------------------

cdef class SeqAlignScore(Serial):

    @staticmethod
    cdef SeqAlignScore _wrap(CRef[CScore] ref):
        cdef SeqAlignScore score = SeqAlignScore.__new__(SeqAlignScore)
        score._ref = ref
        return score

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __iter__(self):
        yield self.id
        yield self.value

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(id={self.id!r}, value={self.value!r})"

    @property
    def id(self):
        if not self._ref.GetNonNullPointer().IsSetId():
            return None
        id_ = &self._ref.GetNonNullPointer().GetIdMut()
        cref = CRef[CObject_id](id_)
        return ObjectId._wrap(cref)

    @property
    def value(self):
        value = &self._ref.GetNonNullPointer().GetValueMut()
        kind = value.Which()
        if kind == CScore_value_choice.e_Int:
            return value.GetInt()
        elif kind == CScore_value_choice.e_Real:
            return value.GetReal()
        raise TypeError(f"Unknown value type: {kind!r}")


cdef class AlignRow:

    @property
    def start(self):
        cdef CSeq_align* obj = self._ref.GetNonNullPointer()
        return obj.GetSeqStart(self._row)

    @property
    def stop(self):
        cdef CSeq_align* obj = self._ref.GetNonNullPointer()
        return obj.GetSeqStop(self._row)

    @property
    def id(self):
        cdef CSeq_align*    obj = self._ref.GetNonNullPointer()
        cdef const CSeq_id* id_ = &obj.GetSeq_id(self._row)
        return SeqId._wrap(CRef[CSeq_id](<CSeq_id*> id_))

    @property
    def num_gap_openings(self):
        cdef CSeq_align* obj = self._ref.GetNonNullPointer()
        return obj.GetNumGapOpenings(self._row)

    @property
    def total_gap_count(self):
        cdef CSeq_align* obj = self._ref.GetNonNullPointer()
        return obj.GetTotalGapCount(self._row)


# --- AlignSegments ------------------------------------------------------------

cdef class AlignSegments(Serial):

    @staticmethod
    cdef AlignSegments _wrap(CRef[C_Segs] ref):
        cdef AlignSegments obj
        cdef C_Segs_choice kind = ref.GetNonNullPointer().Which()

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
        cdef CDense_seg* seg = &self._ref.GetNonNullPointer().GetDensegMut()
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
        cdef CDense_seg* obj = self._ref.GetNonNullPointer()
        return obj.GetNumseg()

    @property
    def ids(self):
        """`list` of `SeqId`: The identifiers of the sequences in the segment.
        """
        cdef CDense_seg* obj = self._ref.GetNonNullPointer()
        cdef list        ids = []

        for ref in obj.GetIdsMut(): # FIXME: const iteration
            ids.append(SeqId._wrap(ref))
        return ids

    @property
    def lengths(self):
        """`list` of `int`: The lengths of each segment.
        """
        # FIXME: make zero-copy
        cdef CDense_seg* obj  = self._ref.GetNonNullPointer()
        cdef list        lens = []

        for length in obj.GetLensMut(): # FIXME: const iteration
            lens.append(length)
        return lens

    @property
    def starts(self):
        """`list` of `int`: The start offsets for the sequences in each segment.
        """
        # FIXME: make zero-copy
        cdef CDense_seg* obj    = self._ref.GetNonNullPointer()
        cdef list        starts = []

        for start in obj.GetStartsMut(): # FIXME: const iteration
            starts.append(start)
        return starts

    @property
    def strands(self):
        """`list` of `str`, or `None`: The strand for each sequence, if any.
        """
        cdef ENa_strand  strand
        cdef CDense_seg* obj     = self._ref.GetNonNullPointer()
        cdef list        strands = []

        for strand in obj.GetStrandsMut():  # FIXME: const iteration
            strands.append(_NA_STRAND_STR[strand])
        return strands

# --- SeqAlign -----------------------------------------------------------------

cdef class SeqAlign(Serial):
    """A sequence alignment, mapping the coordinates of a `BioSeq` to others.

    Sequence alignments are composed of segments, i.e aligned regions which
    contain one or more sequences.

    """

    # FIXME: Reorganize and maybe make the default `__getitem__` access the
    #        segments rather than the rows

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
    def bitscore(self):
        """`float` or `None`: The BLAST-specific bit score.
        """
        cdef double bitscore = NAN
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_BitScore, bitscore):
            return None
        return bitscore

    @property
    def percent_identity(self):
        r"""`float` or `None`: The percent identity of the alignment.

        This refers to the BLAST-style identity, which is computed
        including gaps over the whole alignment:

        .. math::

            id = 100 \times \frac{\text{matches}}{\text{alignment length}}

        """
        cdef int          nident
        cdef unsigned int length = 0
        cdef double       value  = 0

        # return if the score is already available
        if self._ref.GetObject().GetNamedScore(EScoreType.eScore_PercentIdentity, value):
            return value

        # can't compute percent identity if identity is unknown by now
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_IdentityCount, nident):
            return None

        # NOTE: compute the length from the segments, otherwise the alignment length
        # seeems to be inaccurate (see `CAlignFormatUtil::GetPercentIdentity`)
        for l in self._ref.GetObject().GetSegsMut().GetDensegMut().GetLensMut():
            length += l

        if length > 0:
            value = 100.0 * (<double> nident) / (<double> length)

        self._ref.GetObject().SetNamedScore(EScoreType.eScore_PercentIdentity, value)
        return value

    @property
    def percent_coverage(self):
        r"""`float` or `None`: The percent query coverage, if any.

        BLAST ignores the polyA tail when computing coverage for nucleotide
        sequences:

        .. math::

            cov = 100 \times \frac{\text{matches} + \text{mismatches}}
                                  {\text{alignment length} - \text{polyA tail}}

        """
        cdef double cov = NAN
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_PercentCoverage, cov):
            return None
        return cov

    @property
    def evalue(self):
        """`float` or `None`: The BLAST-specific expectation value.
        """
        cdef double evalue = NAN
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_EValue, evalue):
            return None
        return evalue

    @property
    def scores(self):
        # FIXME: Turn into a dict?
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
        cdef CSeq_align*  obj = self._ref.GetNonNullPointer()
        cdef CRef[C_Segs] ref = CRef[C_Segs](&obj.GetSegsMut())
        return AlignSegments._wrap(ref)

    @property
    def alignment_length(self):
        """`int`: The gapped alignment length.
        """
        cdef size_t       length = 0
        cdef CSeq_align*  obj = self._ref.GetNonNullPointer()
        # cdef int alignment_length = obj.GetAlignLength()
        # NOTE: compute the length from the segments, otherwise the alignment length
        # seeems to be inaccurate (see `CAlignFormatUtil::GetPercentIdentity`)
        for l in obj.GetSegsMut().GetDensegMut().GetLensMut():
            length += l
        return length

        return sum(self.segments.data.lengths)

    @property
    def matches(self):
        cdef int nident
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_IdentityCount, nident):
            return None
        return nident

    @property
    def mismatches(self):
        cdef int mm
        if not self._ref.GetObject().GetNamedScore(EScoreType.eScore_MismatchCount, mm):
            return None
        return mm

    @property
    def num_gap_openings(self):
        cdef CSeq_align* obj = self._ref.GetNonNullPointer()
        return obj.GetNumGapOpenings(-1)

    @property
    def total_gap_count(self):
        cdef CSeq_align* obj = self._ref.GetNonNullPointer()
        return obj.GetTotalGapCount(-1)


cdef class GlobalSeqAlign(SeqAlign):
    """A global alignment over the complete lengths of several `BioSeq`.
    """

cdef class DiagonalSeqAlign(SeqAlign):
    pass

cdef class PartialSeqAlign(SeqAlign):
    pass

cdef class DiscontinuousSeqAlign(SeqAlign):
    pass

cdef class SeqAlignSet(Serial):
    """A set of sequence alignments.
    """

    @staticmethod
    cdef SeqAlignSet _wrap(CRef[CSeq_align_set] ref):
        cdef SeqAlignSet obj = SeqAlignSet.__new__(SeqAlignSet)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __iter__(self):
        cdef CRef[CSeq_align] ref
        for ref in self._ref.GetNonNullPointer().GetMut():
            yield SeqAlign._wrap(ref)

    def __len__(self):
        return self._ref.GetNonNullPointer().Get().size()

    cpdef int master_coverage(self) except? 0:
        """Compute the master coverage of the alignment set.
        """
        return CAlignFormatUtil.GetMasterCoverage(self._ref.GetObject())
