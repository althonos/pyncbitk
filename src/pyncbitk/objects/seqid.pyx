# cython: language_level=3

from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.general.dbtag cimport CDbtag
from ..toolkit.objects.general.object_id cimport CObject_id
from ..toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc, E_Choice as CSeq_loc_choice
from ..toolkit.objects.seqloc.seq_interval cimport CSeq_interval
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id, E_Choice as CSeq_id_choice, E_SIC as CSeq_id_SIC
from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.ncbimisc cimport TSeqPos
from ..toolkit.corelib.tempstr cimport CTempString

from ..serial cimport Serial
from .general cimport ObjectId, DBTag

# --- SeqId --------------------------------------------------------------------

cdef class SeqId(Serial):
    """An abstract base class for defining a sequence identifier.
    """

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
        elif kind == CSeq_id_choice.e_Other:
            obj = OtherId.__new__(OtherId)
        else:
            raise NotImplementedError(f"{kind!r}")

        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __eq__(self, object other):
        cdef SeqId       other_
        cdef CSeq_id_SIC result

        if not isinstance(other, SeqId):
            return NotImplemented

        other_ = other
        result = self._ref.GetNonNullPointer().Compare(other_._ref.GetObject())
        
        if result == CSeq_id_SIC.e_DIFF:
            return NotImplemented
        elif result == CSeq_id_SIC.e_error:
            raise RuntimeError(f"Failed to compare {self!r} and {other!r}")
        elif result == CSeq_id_SIC.e_YES:
            return True
        else:
            return False

    def __hash__(self):
        return super().__hash__()

    @staticmethod
    def parse(str text not None):
        """Parse an identifier from an arbitrary string.

        Returns:
            `SeqId`: The appropriate `SeqId` subclass for the given 
            identifier string.

        Example:
            >>> SeqId.parse("JH379476.1")
            GenBankId(TextSeqId('JH379476', version=1))

        """
        cdef bytes _text = text.encode()
        cdef CSeq_id* _id = new CSeq_id(CTempString(_text))
        return SeqId._wrap(CRef[CSeq_id](_id))

cdef class LocalId(SeqId):
    """A local identifier for naming privately maintained data.
    """

    def __init__(self, ObjectId id not None):
        """__init__(self, id)\n--\n

        Create a new local identifier.

        Arguments:
            id (`~pyncbitk.objects.general.ObjectId`): The object identifier.

        """
        cdef CSeq_id* obj = new CSeq_id()
        obj.Select(CSeq_id_choice.e_Local)
        obj.SetLocal(id._ref.GetObject())
        self._ref.Reset(obj)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.object_id!r})"

    def __rich_repr__(self):
        yield self.object_id

    def __str__(self):
        return str(self.object_id)

    def __reduce__(self):
        return type(self), (self.object_id,)

    @property
    def object_id(self):
        """`~pyncbitk.objects.general.ObjectId`: The object identifier.
        """
        cdef CObject_id* id = &self._ref.GetNonNullPointer().GetLocalMut()
        return ObjectId._wrap(CRef[CObject_id](id))

cdef class RefSeqId(SeqId):
    """A sequence identifier from the NCBI Reference Sequence project.
    """

cdef class GenBankId(SeqId):
    """A sequence identifier from the NCBI GenBank database.
    """

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    def __rich_repr__(self):
        yield self.id

    @property
    def id(self):
        """`~pyncbitk.objects.seqid.TextSeqId`: The text identifier.
        """
        cdef CTextseq_id* id = &self._ref.GetNonNullPointer().GetGenbankMut()
        return TextSeqId._wrap(CRef[CTextseq_id](id))


cdef class ProteinDataBankId(SeqId):
    """A sequence identifier from the Protein Data Bank.
    """

cdef class GeneralId(SeqId):
    """A sequence identifier from a database.
    """

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    def __rich_repr__(self):
        yield self.db_tag

    @property
    def db_tag(self):
        """`~pyncbitk.objects.general.DBTag`: A database tag.
        """
        cdef CDbtag* id = &self._ref.GetNonNullPointer().GetGeneralMut()
        return DBTag._wrap(CRef[CDbtag](id))


cdef class OtherId(SeqId):
    """A sequence identifier for other databases.
    """

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id!r})"

    def __rich_repr__(self):
        yield self.id

    @property
    def id(self):
        """`~pyncbitk.objects.seqid.TextSeqId`: The text identifier.
        """
        cdef CTextseq_id* id = &self._ref.GetNonNullPointer().GetOtherMut()
        return TextSeqId._wrap(CRef[CTextseq_id](id))

# --- TextSeqId ----------------------------------------------------------------

cdef class TextSeqId(Serial):
    """A text identifier.

    This format is a standardized identifier for different databases (such 
    as GenBank), where the general format is ``<accession>.<version>``.

    """
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

    def __rich_repr__(self):
        yield self.accession
        yield "name", self.name, None
        yield "version", self.version, 0
        yield "release", self.release, None

    @property
    def accession(self):
        """`str` or `None`: The identifier accession, if any.
        """
        if not self._ref.GetNonNullPointer().IsSetAccession():
            return None
        return self._ref.GetNonNullPointer().GetAccession().decode()

    @property
    def name(self):
        """`str` or `None`: The identifier name, if any.
        """
        if not self._ref.GetNonNullPointer().IsSetName():
            return None
        return self._ref.GetNonNullPointer().GetName().decode()

    @property
    def version(self):
        """`str` or `None`: The identifier version, if any.
        """
        if not self._ref.GetNonNullPointer().IsSetVersion():
            return None
        return self._ref.GetNonNullPointer().GetVersion()

    @version.setter
    def version(self, int version):
        self._ref.GetNonNullPointer().SetVersion(version)

    @property
    def release(self):
        """`str` or `None`: The identifier release, if any.
        """
        if not self._ref.GetNonNullPointer().IsSetRelease():
            return None
        return self._ref.GetNonNullPointer().GetRelease().decode()
