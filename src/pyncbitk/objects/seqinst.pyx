# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.seq_inst cimport CSeq_inst, ETopology, EStrand
from ..toolkit.objects.seq.seq_inst cimport EMol as CSeq_inst_mol
from ..toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from ..toolkit.objects.seq.seq_data cimport CSeq_data
from ..toolkit.serial.serialdef cimport ESerialRecursionMode, ESerialDataFormat, ESerial_Xml_Flags

from ..serial cimport Serial
from .seqdata cimport SeqData, SeqNaData, SeqAaData

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

cdef dict _SEQINST_STRAND_STR = {
    EStrand.eStrand_not_set: None,
    EStrand.eStrand_ss: "single",
    EStrand.eStrand_ds: "double",
    EStrand.eStrand_mixed: "mixed",
    EStrand.eStrand_other: "other",
}

cdef dict _SEQINST_STRAND_ENUM = {
    v:k for k,v in _SEQINST_STRAND_STR.items()
}

cdef class SeqInst(Serial):

    @staticmethod
    cdef SeqInst _wrap(CRef[CSeq_inst] ref):
        cdef SeqInst        obj
        cdef CSeq_inst_repr kind = ref.GetObject().GetRepr()

        if kind == CSeq_inst_repr.eRepr_not_set:
            obj = EmptyInst.__new__(EmptyInst)
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
        str strand=None,
        str molecule=None,
        object length=None,
    ):
        # try to detect molecule if possible
        if molecule not in _SEQINST_MOLECULE_ENUM:
            raise ValueError(f"invalid molecule: {molecule!r}")

        # try to detect strand if possible
        if strand is None:
            if molecule == "dna":
                strand = "double"
            elif molecule == "rna" or molecule == "protein":
                strand = "single"
        elif strand not in _SEQINST_STRAND_ENUM:
            raise ValueError(f"invalid strand: {strand!r}")

        # check topology
        if topology not in _SEQINST_TOPOLOGY_ENUM:
            raise ValueError(f"invalid topology: {topology!r}")

        # set data
        cdef CSeq_inst* obj = new CSeq_inst()
        obj.SetMol(_SEQINST_MOLECULE_ENUM[molecule])
        obj.SetStrand(_SEQINST_STRAND_ENUM[strand])
        obj.SetTopology(_SEQINST_TOPOLOGY_ENUM[topology])
        if length is not None:
            obj.SetLength(length)
        self._ref.Reset(obj)

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = []

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strand is not None:
            args.append(f"strand={self.strand!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")

        return f"{ty}({', '.join(args)})"

    @property
    def length(self):
        if not self._ref.GetObject().IsSetLength():
            return None
        return self._ref.GetObject().GetLength()

    @property
    def molecule(self):
        cdef CSeq_inst_mol kind = self._ref.GetPointer().GetMol()
        return _SEQINST_MOLECULE_STR[kind]

    @property
    def topology(self):
        if not self._ref.GetObject().IsSetTopology():
            return None
        return _SEQINST_TOPOLOGY_STR[self._ref.GetObject().GetTopology()]

    @property
    def strand(self):
        if not self._ref.GetObject().IsSetStrand():
            return None
        return _SEQINST_STRAND_STR[self._ref.GetObject().GetStrand()]

    @property
    def data(self):
        if not self._ref.GetObject().IsSetSeq_data():
            return None
        return SeqData._wrap(CRef[CSeq_data](&self._ref.GetObject().GetSeq_dataMut()))


cdef class EmptyInst(SeqInst):
    pass

cdef class VirtualInst(SeqInst):
    pass

cdef class ContinuousInst(SeqInst):

    def __init__(
        self,
        SeqData data,
        *,
        topology="linear",
        strand=None,
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
            strand=strand,
            molecule=molecule,
            length=length
        )

        cdef CSeq_inst* obj = &self._ref.GetObject()
        obj.SetRepr(CSeq_inst_repr.eRepr_raw)
        obj.SetSeq_data(data._ref.GetObject())

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = [repr(self.data)]

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strand is not None:
            args.append(f"strand={self.strand!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")

        return f"{ty}({', '.join(args)})"

cdef class SegmentedInst(SeqInst):
    pass

cdef class ConstructedInst(SeqInst):
    pass

cdef class RefInst(SeqInst):
    pass

cdef class ConsensusInst(SeqInst):
    pass

cdef class MapInst(SeqInst):
    pass

cdef class DeltaInst(SeqInst):

    cpdef ContinuousInst to_continuous(self):
        cdef CSeq_inst* copy = new CSeq_inst()
        copy.Assign(self._ref.GetNonNullPointer()[0], ESerialRecursionMode.eRecursive)
        if not copy.ConvertDeltaToRaw():
            raise ValueError("Could not convert delta instance to continuous")
        return SeqInst._wrap(CRef[CSeq_inst](copy))
