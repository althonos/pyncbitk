# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_inst cimport CSeq_inst, ETopology, EStrand
from ..toolkit.objects.seq.seq_inst cimport EMol as CSeq_inst_mol
from ..toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from ..toolkit.objects.seq.seq_data cimport CSeq_data
from ..toolkit.serial.serialdef cimport ESerialRecursionMode, ESerialDataFormat, ESerial_Xml_Flags

from ..serial cimport Serial
from .seqdata cimport SeqData, SeqNaData, SeqAaData

from .seqdata import SeqData


# --- SeqInst ------------------------------------------------------------------

cdef dict _SEQINST_MOLECULE_STR = {
    CSeq_inst_mol.eMol_not_set: None,
    CSeq_inst_mol.eMol_dna: "dna",
    CSeq_inst_mol.eMol_rna: "rna",
    CSeq_inst_mol.eMol_aa: "protein",
    CSeq_inst_mol.eMol_na: "nucleotide",
    CSeq_inst_mol.eMol_other: "other",
}

cdef dict _SEQINST_MOLECULE_ENUM = {
    v:k for k,v in _SEQINST_MOLECULE_STR.items()
}

cdef dict _SEQINST_TOPOLOGY_STR = {
    ETopology.eTopology_not_set: None,
    ETopology.eTopology_linear: "linear",
    ETopology.eTopology_circular: "circular",
    ETopology.eTopology_tandem: "tandem",
    ETopology.eTopology_other: "other",
}

cdef dict _SEQINST_TOPOLOGY_ENUM = {
    v:k for k,v in _SEQINST_TOPOLOGY_STR.items()
}

cdef dict _SEQINST_STRANDEDNESS_STR = {
    EStrand.eStrand_not_set: None,
    EStrand.eStrand_ss: "single",
    EStrand.eStrand_ds: "double",
    EStrand.eStrand_mixed: "mixed",
    EStrand.eStrand_other: "other",
}

cdef dict _SEQINST_STRANDEDNESS_ENUM = {
    v:k for k,v in _SEQINST_STRANDEDNESS_STR.items()
}

cdef class SeqInst(Serial):
    """Abstract base class for declaring the contents of a sequence.
    """

    @staticmethod
    cdef SeqInst _wrap(CRef[CSeq_inst] ref):
        cdef SeqInst        obj
        cdef CSeq_inst_repr kind = ref.GetNonNullPointer().GetRepr()

        if kind == CSeq_inst_repr.eRepr_not_set:
            obj = SeqInst.__new__(SeqInst)
        elif kind == CSeq_inst_repr.eRepr_virtual:
            obj = VirtualInst.__new__(VirtualInst)
        elif kind == CSeq_inst_repr.eRepr_raw:
            obj = ContinuousInst.__new__(ContinuousInst)
        elif kind == CSeq_inst_repr.eRepr_seg:
            obj = SegmentedInst.__new__(SegmentedInst)
        elif kind == CSeq_inst_repr.eRepr_const:
            obj = ConstructedInst.__new__(ConstructedInst)
        elif kind == CSeq_inst_repr.eRepr_ref:
            obj = RefInst.__new__(RefInst)
        elif kind == CSeq_inst_repr.eRepr_consen:
            obj = ConsensusInst.__new__(ConsensusInst)
        elif kind == CSeq_inst_repr.eRepr_map:
            obj = MapInst.__new__(MapInst)
        elif kind == CSeq_inst_repr.eRepr_delta:
            obj = DeltaInst.__new__(DeltaInst)
        else:
            raise NotImplementedError

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(
        self,
        *,
        str topology="linear",
        str strandedness=None,
        str molecule=None,
        object length=None,
    ):
        # try to detect molecule if possible
        if molecule not in _SEQINST_MOLECULE_ENUM:
            raise ValueError(f"invalid molecule: {molecule!r}")

        # try to detect strandedness if possible
        if strandedness is None:
            if molecule == "dna":
                strandedness = "double"
            elif molecule == "rna" or molecule == "protein":
                strandedness = "single"
        elif strandedness not in _SEQINST_STRANDEDNESS_ENUM:
            raise ValueError(f"invalid strandedness: {strandedness!r}")

        # check topology
        if topology not in _SEQINST_TOPOLOGY_ENUM:
            raise ValueError(f"invalid topology: {topology!r}")

        # set data
        cdef CSeq_inst* obj = new CSeq_inst()
        obj.SetRepr(CSeq_inst_repr.eRepr_not_set)
        obj.SetMol(_SEQINST_MOLECULE_ENUM[molecule])
        obj.SetStrand(_SEQINST_STRANDEDNESS_ENUM[strandedness])
        obj.SetTopology(_SEQINST_TOPOLOGY_ENUM[topology])
        if length is not None:
            obj.SetLength(length)
        self._ref.Reset(obj)

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = []

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strandedness is not None:
            args.append(f"strand={self.strandedness!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")

        return f"{ty}({', '.join(args)})"

    @property
    def length(self):
        """`int`: The length of the sequence.
        """
        if not self._ref.GetObject().IsSetLength():
            return None
        return self._ref.GetObject().GetLength()

    @property
    def molecule(self):
        """`str` or `None`: The kind of molecule of the sequence, if any.
        """
        cdef CSeq_inst_mol kind = self._ref.GetPointer().GetMol()
        return _SEQINST_MOLECULE_STR[kind]

    @property
    def topology(self):
        """`str` or `None`: The topology of the sequence, if any.
        """
        if not self._ref.GetObject().IsSetTopology():
            return None
        return _SEQINST_TOPOLOGY_STR[self._ref.GetNonNullPointer().GetTopology()]

    @property
    def strandedness(self):
        """`str` or `None`: The strandedness of the sequence, if any.
        """
        if not self._ref.GetObject().IsSetStrand():
            return None
        return _SEQINST_STRANDEDNESS_STR[self._ref.GetNonNullPointer().GetStrand()]

    @property
    def data(self):
        """`SeqData` or `None`: The concrete sequence data.
        """
        if not self._ref.GetObject().IsSetSeq_data():
            return None
        return SeqData._wrap(CRef[CSeq_data](&self._ref.GetNonNullPointer().GetSeq_dataMut()))


cdef class VirtualInst(SeqInst):
    """An instance corresponding to a sequence with no data.
    """

cdef class ContinuousInst(SeqInst):
    """An instance corresponding to a single continuous sequence.
    """

    def __init__(
        self,
        SeqData data,
        *,
        topology="linear",
        strandedness=None,
        molecule=None,
        length=None,
    ):
        if length is None and hasattr(data, "length"):
            length = data.length

        if molecule is None:
            if isinstance(data, SeqNaData):
                molecule = "dna"
            elif isinstance(data, SeqAaData):
                molecule = "aa"

        super().__init__(
            topology=topology,
            strandedness=strandedness,
            molecule=molecule,
            length=length
        )

        cdef CSeq_inst* obj = self._ref.GetNonNullPointer()
        obj.SetRepr(CSeq_inst_repr.eRepr_raw)
        obj.SetSeq_data(data._ref.GetObject())

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = [repr(self.data)]

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strandedness is not None:
            args.append(f"strandedness={self.strandedness!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")

        return f"{ty}({', '.join(args)})"

cdef class SegmentedInst(SeqInst):
    """An instance corresponding to a segmented sequence.
    """

cdef class ConstructedInst(SeqInst):
    """An instance corresponding to a constructed sequence.
    """

cdef class RefInst(SeqInst):
    """An instance corresponding to a reference to another sequence.
    """

cdef class ConsensusInst(SeqInst):
    """An instance corresponding to a consensus sequence.
    """

cdef class MapInst(SeqInst):
    """An instance corresponding to an ordered mapping.
    """

cdef class DeltaInst(SeqInst):
    """An instance corresponding to changed applied to other sequences.
    """

    cpdef ContinuousInst to_continuous(self):
        """Transform this instance to a continuous sequence instance.

        Returns:
            `ContinuousInst`: The equivalent instance as a single
            continuous sequence instance.

        """
        cdef CSeq_inst* copy = new CSeq_inst()
        copy.Assign(self._ref.GetNonNullPointer()[0], ESerialRecursionMode.eRecursive)
        if not copy.ConvertDeltaToRaw():
            raise ValueError("Could not convert delta instance to continuous")
        return SeqInst._wrap(CRef[CSeq_inst](copy))
