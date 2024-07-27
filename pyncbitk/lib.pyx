# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.utility cimport move

from .toolkit.algo.blast.api.bl2seq cimport CBl2Seq
from .toolkit.algo.blast.api.blast_types cimport EProgram, ProgramNameToEnum, TSeqAlignVector
from .toolkit.algo.blast.api.sseqloc cimport SSeqLoc, TSeqLocVector
from .toolkit.corelib.ncbiobj cimport CConstRef, CRef
from .toolkit.corelib.ncbistr cimport kEmptyStr
from .toolkit.corelib.ncbitype cimport Uint4
from .toolkit.corelib.tempstr cimport CTempString
from .toolkit.objects.blastdb.blast_def_line_set cimport CBlast_def_line_set
from .toolkit.objects.general.object_id cimport CObject_id
from .toolkit.objects.general.object_id cimport E_Choice as CObject_id_choice
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objects.seq.bioseq cimport ELabelType as CBioseq_ELabelType
from .toolkit.objects.seq.seq_data cimport CSeq_data
from .toolkit.objects.seq.seq_data cimport E_Choice as CSeq_data_choice
from .toolkit.objects.seq.seq_inst cimport CSeq_inst, ETopology, EStrand
from .toolkit.objects.seq.seq_inst cimport EMol as CSeq_inst_mol
from .toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from .toolkit.objects.seq.iupacna cimport CIUPACna
from .toolkit.objects.seqloc.seq_id cimport CSeq_id
from .toolkit.objects.seqloc.seq_id cimport E_Choice as CSeq_id_choice
from .toolkit.objects.seqloc.seq_loc cimport CSeq_loc
from .toolkit.objects.seqloc.seq_loc cimport E_Choice as CSeq_loc_choice
from .toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from .toolkit.objects.seqset.seq_entry cimport CSeq_entry
from .toolkit.objects.seqset.seq_entry cimport E_Choice as CSeq_entry_choice
from .toolkit.objects.seqalign.seq_align cimport CSeq_align
from .toolkit.objects.seqalign.seq_align_set cimport CSeq_align_set
from .toolkit.objects.seqalign.score cimport CScore, C_Value as CScore_value, E_Choice as CScore_value_choice
from .toolkit.objmgr.object_manager cimport CObjectManager
from .toolkit.objmgr.scope cimport CScope
from .toolkit.objtools.readers.fasta cimport CFastaReader
from .toolkit.serial.serialbase cimport CSerialObject, MSerial_Format_AsnText
from .toolkit.serial.serialdef cimport ESerialRecursionMode
from .toolkit.objects.blastdb.blast_def_line cimport CBlast_def_line
from .toolkit.objects.blastdb.blast_def_line_set cimport CBlast_def_line_set

import os

# --- ObjectManager ------------------------------------------------------------

cdef class ObjectManager:
    cdef CRef[CObjectManager] _mgr

    def __init__(self):
        self._mgr = CObjectManager.GetInstance()

    cpdef Scope scope(self):
        return Scope(self)


cdef class Scope:
    cdef CRef[CScope] _scope

    def __init__(self, ObjectManager manager):
        self._scope.Reset(new CScope(manager._mgr.GetObject()))

    def __enter__(self):
        return self

    def __exit__(self, exc_ty, exc_val, traceback):
        self.close()
        return False

    def close(self):
        self._scope.ReleaseOrNull()

    cpdef void add_bioseq(self, BioSeq seq) except *:
        if self._scope.Empty():
            raise RuntimeError("attempted to use a closed scope")

        self._scope.GetPointer().AddBioseq(seq._ref.GetObject())



cdef CBlast_def_line_set* _x = new CBlast_def_line_set()
    


# --- BLAST --------------------------------------------------------------------

cdef class BlastSeqLoc:
    cdef SSeqLoc _seqloc

    def __init__(self, SeqLoc loc, Scope scope, SeqLoc mask = None):
        self._seqloc = SSeqLoc(loc._loc.GetObject(), scope._scope.GetObject())


cdef class Blast:
    cdef CRef[CBl2Seq] _blast

    def __init__(
        self, 
        BlastSeqLoc query, 
        object subject, 
        str program
    ): 

        cdef EProgram      p
        cdef BlastSeqLoc   s
        cdef TSeqLocVector subjects = TSeqLocVector()
        cdef bytes         _program = program.encode()

        try:
            p = ProgramNameToEnum(_program)
        except Exception as e:
            raise ValueError(f"Invalid BLAST program: {program!r}")

        if isinstance(subject, BlastSeqLoc):
            subject = (subject, )

        for s in subject:
            subjects.push_back(s._seqloc)

        


        self._blast.Reset(
            new CBl2Seq(query._seqloc, subjects, p)
        )

    def run(self):
        cdef TSeqAlignVector _raw_alignments 
        cdef SeqAlignSet     ali
        cdef list            alignments      
        
        with nogil:
            _raw_alignments = self._blast.GetObject().Run()

        alignments = []
        for ref in _raw_alignments:
            alignments.append(SeqAlignSet._wrap(ref))

        return alignments


cdef class SeqAlignScore:
    # TODO: inherit
    cdef CRef[CScore] _ref

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
        id_ = &self._ref.GetObject().GetIdRw()
        cref = CRef[CObject_id](id_)
        return ObjectId._wrap(cref)

    @property
    def value(self):
        value = &self._ref.GetObject().GetValueRw()
        kind = value.Which()
        if kind == CScore_value_choice.e_Int:
            return value.GetInt()
        elif kind == CScore_value_choice.e_Real:
            return value.GetReal()
        raise TypeError(f"Unknown value type: {kind!r}")


cdef class SeqAlign:
    cdef CRef[CSeq_align] _ref

    @staticmethod
    cdef SeqAlign _wrap(CRef[CSeq_align] ref):
        cdef SeqAlign obj = SeqAlign.__new__(SeqAlign)
        obj._ref = ref
        return obj

    def dumps(self):
        cdef string s = string()
        s << <CSerialObject&> self._ref.GetNonNullPointer()[0]
        return s.decode()

    @property
    def scores(self):
        cdef CRef[CScore]  ref
        cdef SeqAlignScore score 
        cdef list          scores = []

        if not self._ref.GetObject().IsSetScore():
            return None

        for ref in self._ref.GetObject().GetScoreRw():
            scores.append(SeqAlignScore._wrap(ref))

        return scores



cdef class SeqAlignSet:
    cdef CRef[CSeq_align_set] _ref

    @staticmethod
    cdef SeqAlignSet _wrap(CRef[CSeq_align_set] ref):
        cdef SeqAlignSet obj = SeqAlignSet.__new__(SeqAlignSet)
        obj._ref = ref
        return obj

    def __iter__(self):
        cdef CRef[CSeq_align] ref
        for ref in self._ref.GetObject().GetRw():
            yield SeqAlign._wrap(ref)

    def __len__(self):
        return self._ref.GetObject().Get().size()

# --- ObjectId -----------------------------------------------------------------

cdef class ObjectId:
    """A basic identifier for any NCBI Toolkit object.
    """
    cdef CRef[CObject_id] _ref

    @staticmethod
    cdef ObjectId _wrap(ref: CRef[CObject_id]):
        cdef ObjectId obj
        cdef CObject_id_choice kind = ref.GetPointer().Which()

        if kind == CObject_id_choice.e_Id:
            obj = IntId.__new__(IntId)
        elif kind == CObject_id_choice.e_Str:
            obj = StrId.__new__(StrId)
        else:
            raise NotImplementedError
        
        obj._ref = ref
        return obj

    def __init__(self):
        raise TypeError("Can't instantiate abstract class ObjectId")

cdef class StrId(ObjectId):
    """An object identifier which is stored as a C++ string.
    """

    def __init__(self, str id):
        cdef bytes _id = id.encode()
        cdef CObject_id* obj = new CObject_id()
        obj.Select(CObject_id_choice.e_Str)
        obj.SetStr(_id)
        self._ref = CRef[CObject_id](obj)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    @property
    def id(self):
        """`str`: The value of the string identifier.
        """
        return self._ref.GetNonNullPointer().GetStr().decode()

cdef class IntId(ObjectId):
    """An object identifier which is stored as an integer.
    """

    def __init__(self, int id):
        cdef CObject_id* obj = new CObject_id()
        obj.Select(CObject_id_choice.e_Id)
        obj.SetId(id)
        self._ref = CRef[CObject_id](obj)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    @property
    def id(self):
        """`int`: The value of the integer identifier.
        """
        return self._ref.GetNonNullPointer().GetId()

# --- TextSeqId ----------------------------------------------------------------

cdef class TextSeqId:
    cdef CRef[CTextseq_id] _ref

    @staticmethod
    cdef TextSeqId _wrap(ref: CRef[CTextseq_id]):
        cdef TextSeqId obj = TextSeqId.__new__(TextSeqId)
        obj._ref = ref
        return obj

    def __init__(
        self,
        str accession,
        *,
        str name = None,
        int version = 0,
        str release = None,
        bool allow_dot_version = True
    ):
        cdef bytes _accession = accession.encode()
        cdef bytes _name = None if name is None else name.encode()
        cdef bytes _release = None if release is None else release.encode()

        # TODO
        self._ref.Reset(new CTextseq_id())
        # self._tid.GetNonNullPointer().Set(
        #     CTempString(_accession),
        # )

    def __repr__(self):
        cdef str ty    = type(self).__name__
        cdef list args = [repr(self.accession)]
        
        name = self.name
        if name is not None:
            args.append(f"name={name!r}")
        version = self.version
        if version != 0:
            args.append(f"version={version!r}")
        release = self.release
        if release is not None:
            args.append(f"release={release!r}")

        return f"{ty}({', '.join(args)})"

    @property
    def accession(self):
        if not self._ref.GetObject().IsSetAccession():
            return None
        return self._ref.GetObject().GetAccession().decode()

    @property
    def name(self):
        if not self._ref.GetObject().IsSetName():
            return None
        return self._ref.GetObject().GetName().decode()

    @property
    def version(self):
        if not self._ref.GetObject().IsSetVersion():
            return None
        return self._ref.GetObject().GetVersion()

    @version.setter
    def version(self, int version):
        self._ref.GetObject().SetVersion(version)

    @property
    def release(self):
        if not self._ref.GetObject().IsSetRelease():
            return None
        return self._ref.GetObject().GetRelease().decode()


# --- SeqId --------------------------------------------------------------------

cdef class SeqId:
    cdef CRef[CSeq_id] _ref

    @staticmethod
    cdef SeqId _wrap(ref: CRef[CSeq_id]):
        cdef SeqId obj
        cdef CSeq_id_choice kind = ref.GetPointer().Which()

        if kind == CSeq_id_choice.e_Local:
            obj = LocalId.__new__(LocalId)
        elif kind == CSeq_id_choice.e_Genbank:
            obj = GenBankId.__new__(GenBankId)
        else:
            raise NotImplementedError
        
        obj._ref = ref
        return obj

    @staticmethod
    def parse(str text):
        cdef bytes _text = text.encode()
        cdef CSeq_id* _id = new CSeq_id(CTempString(_text))
        return SeqId._wrap(CRef[CSeq_id](_id))

    # FIXME: debug
    def dumps(self):
        cdef string s = string()
        s << <CSerialObject&> self._ref.GetNonNullPointer()[0]
        return s.decode()


cdef class LocalId(SeqId):

    def __init__(self, ObjectId id):
        cdef CSeq_id* obj = new CSeq_id()
        obj.Select(CSeq_id_choice.e_Local)
        obj.SetLocal(id._ref.GetObject())
        self._ref.Reset(obj)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"
    
    @property
    def id(self):
        cdef CObject_id* id = &self._ref.GetNonNullPointer().GetLocalRw()
        return ObjectId._wrap(CRef[CObject_id](id))
        
cdef class RefSeqId(SeqId):
    pass

cdef class GenBankId(SeqId):
    
    @property
    def id(self):
        cdef CTextseq_id* id = &self._ref.GetObject().GetGenbankRw()
        return TextSeqId._wrap(CRef[CTextseq_id](id))


cdef class ProteinDataBankId(SeqId):
    pass

cdef class GeneralId(SeqId):
    pass

# --- BioSeq -------------------------------------------------------------------

cdef class BioSeq:
    cdef CRef[CBioseq] _ref

    @staticmethod
    cdef BioSeq _wrap(ref: CRef[CBioseq]):
        cdef BioSeq obj = BioSeq.__new__(BioSeq)
        obj._ref = ref
        return obj
    
    def __init__(
        self,
        object ids,
        SeqInst instance,
    ):
        cdef SeqId    _id
        cdef CBioseq* _bioseq = new CBioseq()

        _bioseq.SetInst(instance._ref.GetObject())
        for _id in ids:
            _bioseq.SetId().push_back(_id._ref)

        self._ref.Reset(new CBioseq())

    @property
    def ids(self):
        assert self._ref.GetNonNullPointer().IsSetId()  # mandatory

        cdef SeqId                  seqid
        cdef list                   ids   = []
        cdef cpplist[CRef[CSeq_id]] _ids  = self._ref.GetObject().GetId()

        for ref in _ids:
            ids.append(SeqId._wrap(ref))

        return ids

    @property
    def instance(self):
        assert self._ref.GetNonNullPointer().IsSetInst()  # mandatory
        return SeqInst._wrap(CRef[CSeq_inst](&self._ref.GetObject().GetInstRw()))

    # FIXME: debug
    def dumps(self):
        cdef string s = string()
        s << <CSerialObject&> self._ref.GetNonNullPointer()[0]
        return s.decode()
        

# --- SeqInst ------------------------------------------------------------------

cdef class SeqInst:
    cdef CRef[CSeq_inst] _ref

    @staticmethod
    cdef SeqInst _wrap(ref: CRef[CSeq_inst]):
        cdef SeqInst obj
        cdef CSeq_inst_repr kind = ref.GetPointer().GetRepr()

        if kind == CSeq_inst_repr.eRepr_not_set:
            obj = EmptyInst.__new__(LocalId)
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

    @property
    def length(self):
        if not self._ref.GetObject().IsSetLength():
            return None
        return self._ref.GetObject().GetLength()

    @property
    def molecule(self):
        cdef CSeq_inst_mol kind = self._ref.GetPointer().GetMol()
        return CSeq_inst.GetMoleculeClass(kind).decode()

    @property
    def topology(self):
        if not self._ref.GetObject().IsSetTopology():
            return None

        cdef ETopology topology = self._ref.GetObject().GetTopology()
        if topology == ETopology.eTopology_not_set:
            return None
        elif topology == ETopology.eTopology_linear:
            return "linear"
        elif topology == ETopology.eTopology_linear:
            return "circular"
        elif topology == ETopology.eTopology_tandem:
            return "tandem"
        elif topology == ETopology.eTopology_other:
            return "other"
        raise ValueError("unexepected topology")

    @property
    def strand(self):
        if not self._ref.GetObject().IsSetStrand():
            return None

        cdef EStrand strand = self._ref.GetObject().GetStrand()
        if strand == EStrand.eStrand_not_set:
            return None
        elif strand == EStrand.eStrand_ss:
            return "single"
        elif strand == EStrand.eStrand_ds:
            return "double"
        elif strand == EStrand.eStrand_mixed:
            return "mixed"
        elif strand == EStrand.eStrand_other:
            return "other"
        raise ValueError("unexepected strand")

    @property
    def data(self):
        if not self._ref.GetObject().IsSetSeq_data():
            return None
        return SeqData._wrap(CRef[CSeq_data](&self._ref.GetObject().GetSeq_dataRw()))


cdef class EmptyInst(SeqInst):
    pass

cdef class VirtualInst(SeqInst):
    pass

cdef class ContinuousInst(SeqInst):
    pass

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


# --- SeqData ------------------------------------------------------------------

cdef class SeqData:
    cdef CRef[CSeq_data] _ref

    @staticmethod
    cdef SeqData _wrap(ref: CRef[CSeq_data]):
        cdef SeqData obj
        cdef CSeq_data_choice kind = ref.GetPointer().Which()

        if kind == CSeq_data_choice.e_Iupacna:
            obj = IupacNaData.__new__(IupacNaData)
        elif kind == CSeq_data_choice.e_Iupacaa:
            obj = IupacAaData.__new__(IupacAaData)
        elif kind == CSeq_data_choice.e_Ncbi2na:
            obj = Ncbi2NaData.__new__(Ncbi2NaData)
        elif kind == CSeq_data_choice.e_Ncbi4na:
            obj = Ncbi4NaData.__new__(Ncbi4NaData)
        elif kind == CSeq_data_choice.e_Ncbi8na:
            obj = Ncbi8NaData.__new__(Ncbi8NaData)
        elif kind == CSeq_data_choice.e_Ncbipna:
            obj = NcbiPNaData.__new__(NcbiPNaData)
        elif kind == CSeq_data_choice.e_Ncbi8aa:
            obj = Ncbi8AaData.__new__(Ncbi8AaData)
        elif kind == CSeq_data_choice.e_Ncbieaa:
            obj = NcbiEAaData.__new__(NcbiEAaData)
        elif kind == CSeq_data_choice.e_Ncbipaa:
            obj = NcbiPAaData.__new__(NcbiPAaData)
        elif kind == CSeq_data_choice.e_Ncbistdaa:
            obj = NcbiStdAa.__new__(NcbiStdAa)
        elif kind == CSeq_data_choice.e_Gap:
            obj = GapData.__new__(GapData)
        else:
            raise NotImplementedError

        obj._ref = ref
        return obj

cdef class IupacNaData(SeqData):

    def __init__(self, object data):
        cdef bytes      _data
        cdef CSeq_data* _raw  = new CSeq_data()

        if isinstance(data, str):
            _data = data.encode()
        else:
            _data = data

        _raw.Select(CSeq_data_choice.e_Iupacna)
        _raw.SetIupacna(CIUPACna(<string> _data))
    
    cpdef str decode(self):
        return self._data.GetNonNullPointer().GetIupacna().Get().decode()


cdef class IupacAaData(SeqData):
    pass

cdef class Ncbi2NaData(SeqData):
    pass

cdef class Ncbi4NaData(SeqData):
    pass

cdef class Ncbi8NaData(SeqData):
    pass

cdef class NcbiPNaData(SeqData):
    pass

cdef class Ncbi8AaData(SeqData):
    pass

cdef class NcbiEAaData(SeqData):
    pass

cdef class NcbiPAaData(SeqData):
    pass

cdef class NcbiStdAa(SeqData):
    pass

cdef class GapData(SeqData):
    pass


# --- Entry --------------------------------------------------------------------

cdef class Entry:
    cdef CRef[CSeq_entry] _ref

    @staticmethod
    cdef Entry _wrap(ref: CRef[CSeq_entry]):
        
        cdef Entry entry
        cdef CSeq_entry_choice kind = ref.GetNonNullPointer().Which()

        if kind == CSeq_entry_choice.e_Seq:
            entry = SeqEntry.__new__(SeqEntry)
        elif kind == CSeq_entry_choice.e_Set:
            entry = SetEntry.__new__(SetEntry)
        else:
            raise RuntimeError("seq entry kind not defined")
        entry._ref = ref
        return entry

cdef class SeqEntry(Entry):

    @property
    def seq(self):
        cdef CBioseq* bioseq = &self._ref.GetNonNullPointer().SetSeq()
        return BioSeq._wrap(CRef[CBioseq](bioseq))


cdef class SetEntry(Entry):
    pass


# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc:
    cdef CRef[CSeq_loc] _loc

    # FIXME: debug
    def dumps(self):
        cdef string s = string()
        s << <CSerialObject&> self._loc.GetNonNullPointer()[0]
        return s.decode()

cdef class NullLoc(SeqLoc):
    pass

cdef class EmptySeqLoc(SeqLoc):
    pass

cdef class WholeSeqLoc(SeqLoc):
    
    def __init__(self, SeqId id):
        cdef CSeq_loc* loc = new CSeq_loc()
        loc.Select(CSeq_loc_choice.e_Whole)
        loc.SetWhole( id._ref.GetNonNullPointer()[0] )
        self._loc.Reset(loc)

    @property
    def id(self):
        id_ = CRef[CSeq_id](&self._loc.GetNonNullPointer().GetWholeRw())
        return SeqId._wrap(id_)

cdef class SeqIntervalLoc(SeqLoc):
    pass

cdef class PackedSeqLoc(SeqLoc):
    pass

cdef class PackedPointsLoc(SeqLoc):
    pass

cdef class MixLoc(SeqLoc):
    pass

cdef class EquivalentLoc(SeqLoc):
    pass

cdef class BondLoc(SeqLoc):
    pass

cdef class FeatureLoc(SeqLoc):
    pass


# --- FastaReader --------------------------------------------------------------

cdef class FastaReader:
    cdef CFastaReader* _reader

    def __cinit__(self):
        self._reader = NULL

    def __dealloc__(self):
        del self._reader

    def __init__(self, path):
        cdef bytes _path = os.fsencode(path)
        self._reader = new CFastaReader(_path)

    def read(self):
        assert self._reader != NULL
        _entry = self._reader.ReadOneSeq()
        return Entry._wrap(_entry)
        