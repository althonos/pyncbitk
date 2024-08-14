# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.general.object_id cimport CObject_id
from ..toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc, E_Choice as CSeq_loc_choice
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id, E_Choice as CSeq_id_choice
from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.tempstr cimport CTempString

from ..serial cimport Serial
from .general cimport ObjectId

# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc(Serial):

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._loc.GetNonNullPointer()

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
        id_ = CRef[CSeq_id](&self._loc.GetNonNullPointer().GetWholeMut())
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

# --- SeqId --------------------------------------------------------------------

cdef class SeqId(Serial):

    @staticmethod
    cdef SeqId _wrap(CRef[CSeq_id] ref):
        cdef SeqId obj
        cdef CSeq_id_choice kind = ref.GetPointer().Which()

        if kind == CSeq_id_choice.e_Local:
            obj = LocalId.__new__(LocalId)
        elif kind == CSeq_id_choice.e_Genbank:
            obj = GenBankId.__new__(GenBankId)
        elif kind == CSeq_id_choice.e_General:
            obj = GeneralId.__new__(GeneralId)
        else:
            raise NotImplementedError(f"{kind!r}")

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    @staticmethod
    def parse(str text):
        cdef bytes _text = text.encode()
        cdef CSeq_id* _id = new CSeq_id(CTempString(_text))
        return SeqId._wrap(CRef[CSeq_id](_id))

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
        cdef CObject_id* id = &self._ref.GetNonNullPointer().GetLocalMut()
        return ObjectId._wrap(CRef[CObject_id](id))

cdef class RefSeqId(SeqId):
    pass

cdef class GenBankId(SeqId):

    @property
    def id(self):
        cdef CTextseq_id* id = &self._ref.GetObject().GetGenbankMut()
        return TextSeqId._wrap(CRef[CTextseq_id](id))


cdef class ProteinDataBankId(SeqId):
    pass

cdef class GeneralId(SeqId):
    pass

cdef class TextSeqId(Serial):
    # FIXME: Consider removing this data class and using instead an abstract
    #        subclass for `SeqId` that exposes the text seq ID attributes for 
    #        the relevant IDs (GenBank ID, etc)?

    @staticmethod
    cdef TextSeqId _wrap(CRef[CTextseq_id] ref):
        cdef TextSeqId obj = TextSeqId.__new__(TextSeqId)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

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
        cdef bytes _name      = None if name is None else name.encode()
        cdef bytes _release   = None if release is None else release.encode()

        # TODO
        self._ref.Reset(new CTextseq_id())
        raise NotImplementedError
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
