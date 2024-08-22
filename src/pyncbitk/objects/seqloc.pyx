# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.general.object_id cimport CObject_id
from ..toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc, E_Choice as CSeq_loc_choice
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id, E_Choice as CSeq_id_choice, E_SIC as CSeq_id_SIC
from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.tempstr cimport CTempString

from ..serial cimport Serial
from .general cimport ObjectId

# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc(Serial):
    """An abstract base class for defining a location in a `BioSeq`.
    """

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._loc.GetNonNullPointer()

cdef class NullLoc(SeqLoc):
    """A gap of unknown size.
    """

cdef class EmptySeqLoc(SeqLoc):
    """A gap of unknown size inside an alignment.
    """

cdef class WholeSeqLoc(SeqLoc):
    """A reference to an entire `BioSeq`.
    """

    def __init__(self, SeqId sequence_id):
        """__init__(self, id)\n--\n

        Create a new location referencing the given sequence.

        Arguments:
            sequence_id (`~pyncbitk.objects.seqloc.SeqId`): The identifier
                of the sequence being referenced. 

        """
        cdef CSeq_loc* loc = new CSeq_loc()
        loc.Select(CSeq_loc_choice.e_Whole)
        loc.SetWhole(sequence_id._ref.GetNonNullPointer()[0] )
        self._loc.Reset(loc)

    @property
    def sequence_id(self):
        """`~pyncbitk.objects.seqloc.SeqId`: The identifier of the sequence.
        """
        id_ = CRef[CSeq_id](&self._loc.GetNonNullPointer().GetWholeMut())
        return SeqId._wrap(id_)

cdef class SeqIntervalLoc(SeqLoc):
    """A reference to an interval on a `BioSeq`.
    """
    
    @property
    def sequence_id(self):
        """`~pyncbitk.objects.seqloc.SeqId`: The identifier of the sequence.
        """
        data = &self._loc.GetNonNullPointer().GetIntMut()
        id_ = CRef[CSeq_id](&data.GetIdMut())
        return SeqId._wrap(id_)

    @property
    def start(self):
        """`int`: The beginining of the sequence interval.
        """
        data = &self._loc.GetNonNullPointer().GetIntMut()
        return data.GetFrom()

    @property
    def stop(self):
        """`int`: The beginining of the sequence interval.
        """
        data = &self._loc.GetNonNullPointer().GetIntMut()
        return data.GetTo()


cdef class PackedSeqLoc(SeqLoc):
    """A reference to a series of intervals on a `BioSeq`.
    """

cdef class PointLoc(SeqLoc):
    """A reference to a single point on a `BioSeq`.
    """

cdef class PackedPointsLoc(SeqLoc):
    """A reference to a series of points on a `BioSeq`.
    """

cdef class MixLoc(SeqLoc):
    """An arbitrarily complex location.
    """

cdef class EquivalentLoc(SeqLoc):
    """A set of equivalent locations.
    """

cdef class BondLoc(SeqLoc):
    """A chemical bond between two residues.
    """

cdef class FeatureLoc(SeqLoc):
    """A location indirectly referenced through a feature.
    """

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

    @property
    def id(self):
        cdef CTextseq_id* id = &self._ref.GetNonNullPointer().GetGenbankMut()
        return TextSeqId._wrap(CRef[CTextseq_id](id))


cdef class ProteinDataBankId(SeqId):
    """A sequence identifier from the Protein Data Bank.
    """

cdef class GeneralId(SeqId):
    """A sequence identifier from a local database.
    """

# --- TextSeqId ----------------------------------------------------------------

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
        if not self._ref.GetNonNullPointer().IsSetAccession():
            return None
        return self._ref.GetNonNullPointer().GetAccession().decode()

    @property
    def name(self):
        if not self._ref.GetNonNullPointer().IsSetName():
            return None
        return self._ref.GetNonNullPointer().GetName().decode()

    @property
    def version(self):
        if not self._ref.GetNonNullPointer().IsSetVersion():
            return None
        return self._ref.GetNonNullPointer().GetVersion()

    @version.setter
    def version(self, int version):
        self._ref.GetNonNullPointer().SetVersion(version)

    @property
    def release(self):
        if not self._ref.GetNonNullPointer().IsSetRelease():
            return None
        return self._ref.GetNonNullPointer().GetRelease().decode()
