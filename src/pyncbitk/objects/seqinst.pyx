# cython: language_level=3
"""Instantiation of biological sequences.

The `BioSeq` model of the NCBI C++ Toolkit requires two components: one
or more sequence identifiers, which are used to name a sequence *instance*.
A sequence instance represents concrete properties of a biological
sequence, such as its molecule type (DNA, RNA, protein) or its strandedness
if applicable.

Note:
    In the original C++ code, sequence instances follow two parallel class
    hierarchies for representation: the molecule type branch, and the
    representation branch. Because Cython does not support multiple
    inheritance, the classes in this module follow only the representation
    branch of the hierarchy. The molecule type, if known, can be accessed
    with the `~SeqInst.molecule` property.

See Also:
    The `Data mode <https://ncbi.github.io/cxx-toolkit/pages/ch_datamod#ch_datamod.Classes_of_Biologica>`_
    chapter of the NCBI C++ Toolkit documentation.

"""

from cython.operator cimport dereference, preincrement
from libcpp.list cimport list as cpplist

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.ncbimisc cimport TSeqPos
from ..toolkit.objects.seq.seq_inst cimport CSeq_inst, ETopology, EStrand
from ..toolkit.objects.seq.seq_inst cimport EMol as CSeq_inst_mol
from ..toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from ..toolkit.objects.seq.seq_ext cimport CSeq_ext, E_Choice as CSeq_ext_choice
from ..toolkit.objects.seq.ref_ext cimport CRef_ext
from ..toolkit.objects.seq.delta_ext cimport CDelta_ext
from ..toolkit.objects.seq.seq_data cimport CSeq_data
from ..toolkit.objects.seq.seq_literal cimport CSeq_literal
from ..toolkit.objects.seq.delta_seq cimport E_Choice as CDelta_seq_choice
from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc
from ..toolkit.serial.serialdef cimport ESerialRecursionMode, ESerialDataFormat, ESerial_Xml_Flags

from ..serial cimport Serial
from .seqloc cimport SeqLoc
from .seqdata cimport SeqData, SeqNaData, SeqAaData

import functools

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
        if self.length is not None:
            args.append(f"length={self.length!r}")

        return f"{ty}({', '.join(args)})"

    def __rich_repr__(self):
        yield "topology", self.topology, "linear"
        yield "strandedness", self.strandedness, None
        yield "molecule", self.molecule, None
        yield "length", self.length, None

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


cdef class VirtualInst(SeqInst):
    """An instance corresponding to a sequence with no data.

    This class allows to describe the properties of a sequence, such as its
    length or its molecule type, without actually knowing the sequence data.

    """

cdef class ContinuousInst(SeqInst):
    """An instance corresponding to a single continuous sequence.

    This class describes the most simple sequence kind, where we actually know
    the sequence data, which can be reached with the `~SeqInst.data`
    property.

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
        """__init__(self, data, *, topology="linear", strandedness=None, molecule=None, length=None)\n--\n

        Create a new continuous instance from the given sequence data.

        Arguments:
            data (`~pyncbitk.objects.seqdata.SeqData`): The concrete sequence
                data.
            topology (`str`): The topology of the sequence; either ``linear``,
                ``circular``, ``tandem``, ``other``, or `None` for unknown
                sequence topologies.
            strandedness (`str`): The strandedness of the sequence; either
                ``single``, ``double``, ``mixed``, ``other``, or `None` for
                unknown strandedness.
            molecule (`str`): The type of molecule described by the data;
                either ``dna``, ``rna``, ``protein``, ``nucleotide``,
                ``other`` or `None` for unknown molecule types.
            length (`int`): The length of the sequence, if known.

        """
        if length is None and hasattr(data, "length"):
            length = data.length

        if molecule is None:
            if isinstance(data, SeqNaData):
                molecule = "dna"
            elif isinstance(data, SeqAaData):
                molecule = "protein"

        super().__init__(
            topology=topology,
            strandedness=strandedness,
            molecule=molecule,
            length=length
        )

        cdef CSeq_inst* obj = self._ref.GetNonNullPointer()
        obj.SetRepr(CSeq_inst_repr.eRepr_raw)
        obj.SetSeq_data(data._ref.GetObject())

    def __reduce__(self):
        return functools.partial(
            type(self),
            self.data,
            topology=self.topology,
            strandedness=self.strandedness,
            molecule=self.molecule,
            length=self.length,
        ), ()

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = [repr(self.data)]

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strandedness is not None:
            args.append(f"strandedness={self.strandedness!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")
        if self.length is not None:
            args.append(f"length={self.length!r}")

        return f"{ty}({', '.join(args)})"

    def __rich_repr__(self):
        yield self.data
        yield from super().__rich_repr__()

    @property
    def data(self):
        """`SeqData` or `None`: The concrete sequence data.
        """
        if not self._ref.GetObject().IsSetSeq_data():
            return None
        return SeqData._wrap(CRef[CSeq_data](&self._ref.GetNonNullPointer().GetSeq_dataMut()))


cdef class SegmentedInst(SeqInst):
    """An instance corresponding to a segmented sequence.
    """

cdef class ConstructedInst(SeqInst):
    """An instance corresponding to a constructed sequence.
    """

cdef class RefInst(SeqInst):
    """An instance corresponding to a reference to another sequence.

    This class allows to describe the sequence data in terms of a location
    in another sequence, described with a `~pyncbitk.objects.seqloc.SeqLoc`.
    This can be used to alias certain regions of a sequence, such as creating
    references to the genes of contig without having to copy the sequence
    data.

    """

    def __init__(
        self,
        SeqLoc seqloc,
        *,
        topology="linear",
        strandedness=None,
        molecule=None,
        length=None,
    ):
        """__init__(self, seqloc, *, topology="linear", strandedness=None, molecule=None, length=None)\n--\n

        Create a new instance referencing the given location.

        Arguments:
            seqloc (`~pyncbitk.objects.seqloc.SeqLoc`): The location of
                the actual sequence data.
            topology (`str`): The topology of the sequence; either ``linear``,
                ``circular``, ``tandem``, ``other``, or `None` for unknown
                sequence topologies.
            strandedness (`str`): The strandedness of the sequence; either
                ``single``, ``double``, ``mixed``, ``other``, or `None` for
                unknown strandedness.
            molecule (`str`): The type of molecule described by the data;
                either ``dna``, ``rna``, ``protein``, ``nucleotide``,
                ``other`` or `None` for unknown molecule types.
            length (`int`): The length of the sequence, if known.

        """
        super().__init__(
            topology=topology,
            strandedness=strandedness,
            molecule=molecule,
            length=length
        )
        # copy the seqloc into the object
        cdef CRef_ext* ref = new CRef_ext()
        ref.Assign(seqloc._loc.GetObject(), ESerialRecursionMode.eRecursive)
        # add the seqloc into the sequence external data
        cdef CSeq_ext* ext = new CSeq_ext()
        ext.SetRef(ref[0])
        # create the sequence instance
        cdef CSeq_inst* obj = self._ref.GetNonNullPointer()
        obj.SetRepr(CSeq_inst_repr.eRepr_ref)
        obj.SetExt(ext[0])

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = [repr(self.seqloc)]

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strandedness is not None:
            args.append(f"strandedness={self.strandedness!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")
        if self.length is not None:
            args.append(f"length={self.length!r}")

        return f"{ty}({', '.join(args)})"

    def __rich_repr__(self):
        yield self.seqloc
        yield from super().__rich_repr__()

    def __reduce__(self):
        return functools.partial(
            type(self),
            self.seqloc,
            topology=self.topology,
            strandedness=self.strandedness,
            molecule=self.molecule,
            length=self.length,
        ), ()

    @property
    def seqloc(self):
        """`~pyncbitk.objects.seqloc.SeqLoc`: The reference sequence location.
        """
        cdef CRef_ext* ext = &self._ref.GetObject().GetExtMut().GetRefMut()
        cdef CRef[CSeq_loc] ref = CRef[CSeq_loc](&ext.GetMut())
        return SeqLoc._wrap(ref)

cdef class ConsensusInst(SeqInst):
    """An instance corresponding to a consensus sequence.
    """

cdef class MapInst(SeqInst):
    """An instance corresponding to an ordered mapping.
    """

cdef class DeltaInst(SeqInst):
    """An instance corresponding to changed applied to other sequences.
    """

    def __init__(
        self,
        object deltas = (),
        *,
        topology="linear",
        strandedness=None,
        molecule=None,
        length=None
    ):
        cdef Delta       delta
        cdef CDelta_ext* delta_ext = new CDelta_ext()

        for delta in deltas:
            delta_ext.GetMut().push_back(delta._ref)

        super().__init__(
            topology=topology,
            strandedness=strandedness,
            molecule=molecule,
            length=length
        )

        # create the sequence extension to store the delta
        cdef CSeq_ext* ext = new CSeq_ext()
        ext.Select(CSeq_ext_choice.e_Delta)
        ext.SetDelta(delta_ext[0])

        # create the sequence instance
        cdef CSeq_inst* obj = self._ref.GetNonNullPointer()
        obj.SetRepr(CSeq_inst_repr.eRepr_delta)
        obj.SetExt(ext[0])

    def __bool__(self):
        return not self._ref.GetObject().GetExt().GetDelta().Get().empty()

    def __len__(self):
        return self._ref.GetObject().GetExt().GetDelta().Get().size()

    def __iter__(self):
        cdef CRef[CDelta_seq]                         delta
        cdef const cpplist[CRef[CDelta_seq]]*         deltas = &self._ref.GetObject().GetExt().GetDelta().Get()
        cdef cpplist[CRef[CDelta_seq]].const_iterator it     = deltas.const_begin()

        while it != deltas.const_end():
            yield Delta._wrap(dereference(it))
            preincrement(it)

    def __reduce__(self):
        return functools.partial(
            type(self),
            list(self),
            topology=self.topology,
            strandedness=self.strandedness,
            molecule=self.molecule,
            length=self.length
        ), ()

    def __repr__(self):
        cdef str ty    = self.__class__.__name__
        cdef list args = [repr([delta for delta in self])]

        if self.topology != "linear":
            args.append(f"topology={self.topology!r}")
        if self.strandedness is not None:
            args.append(f"strandedness={self.strandedness!r}")
        if self.molecule is not None:
            args.append(f"molecule={self.molecule!r}")
        if self.length is not None:
            args.append(f"length={self.length!r}")

        return f"{ty}({', '.join(args)})"

    def __rich_repr__(self):
        yield list(self)
        yield from super().__rich_repr__()

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


cdef class Delta(Serial):
    """A single delta segment in a `DeltaInst` object.
    """

    @staticmethod
    cdef Delta _wrap(CRef[CDelta_seq] ref):
        cdef Delta                    obj
        cdef CDelta_seq_choice kind = ref.GetNonNullPointer().Which()

        if kind == CDelta_seq_choice.e_not_set:
            obj = Delta.__new__(Delta)
        elif kind == CDelta_seq_choice.e_Loc:
            obj = LocDelta.__new__(LocDelta)
        elif kind == CDelta_seq_choice.e_Literal:
            obj = LiteralDelta.__new__(LiteralDelta)
        else:
            raise RuntimeError("Unsupported `Delta` type")

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()


cdef class LiteralDelta(Delta):
    """A literal delta segment.
    """

    # TODO: handle int fuzz

    def __init__(self, TSeqPos length, SeqData data = None):
        cdef CSeq_literal* lit = new CSeq_literal()
        lit.SetLength(length)
        lit.SetSeq_data(data._ref.GetObject())
        cdef CDelta_seq* delta = new CDelta_seq()
        delta.Select(CDelta_seq_choice.e_Literal)
        delta.SetLiteral(lit[0])
        self._ref.Reset(delta)

    def __repr__(self):
        cdef TSeqPos length = self._ref.GetObject().GetLiteral().GetLength()
        cdef SeqData data   = self.data
        cdef list    args   = [repr(length)]
        if data is not None:
            args.append(repr(data))
        return f"{type(self).__name__}({', '.join(args)})"

    def __rich_repr__(self):
        yield self.length
        yield None, self.data, None

    def __reduce__(self):
        return type(self), (self.length, self.data)

    @property
    def length(self):
        """`int`: The length of the literal delta.
        """
        cdef const CSeq_literal* lit = &self._ref.GetObject().GetLiteral()
        return lit.GetLength()

    @property
    def data(self):
        """`SeqData` or `None`: The concrete sequence data, if set.
        """
        cdef const CSeq_literal* lit = &self._ref.GetObject().GetLiteral()
        if not lit.IsSetSeq_data():
            return None
        cdef const CSeq_data* data = &lit.GetSeq_data()
        return SeqData._wrap(CRef[CSeq_data](<CSeq_data*> data))


cdef class LocDelta(Delta):
    """A sequence location delta segment.
    """

    def __init__(self, SeqLoc seqloc not None):
        cdef CDelta_seq* delta = new CDelta_seq()
        delta.Select(CDelta_seq_choice.e_Loc)
        delta.SetLoc(seqloc._loc.GetObject())
        self._ref.Reset(delta)

    def __reduce__(self):
        return type(self), (self.seqloc,)

    def __repr__(self):
        return f"{type(self).__name__}({self.seqloc!r})"

    def __rich_repr__(self):
        yield self.seqloc

    @property
    def seqloc(self):
        """`~pyncbitk.objects.seqloc.SeqLoc`: The location for this delta.
        """
        cdef CRef[CSeq_loc] loc = CRef[CSeq_loc](&self._ref.GetObject().GetLocMut())
        return SeqLoc._wrap(loc)
