# cython: language_level=3
"""Locations on a biological sequence.

This module contains classes that are used to describe a location on a
`BioSeq`. 

"""

from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.general.object_id cimport CObject_id
from ..toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc, E_Choice as CSeq_loc_choice
from ..toolkit.objects.seqloc.seq_interval cimport CSeq_interval
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id, E_Choice as CSeq_id_choice, E_SIC as CSeq_id_SIC
from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.ncbimisc cimport TSeqPos
from ..toolkit.corelib.tempstr cimport CTempString

from ..serial cimport Serial
from .general cimport ObjectId
from .seqid cimport SeqId

# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc(Serial):
    """An abstract base class for defining a location in a `BioSeq`.
    """

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._loc.GetNonNullPointer()

    @staticmethod
    cdef SeqLoc _wrap(CRef[CSeq_loc] ref):
        cdef SeqLoc obj
        cdef CSeq_loc_choice kind = ref.GetPointer().Which()

        if kind == CSeq_loc_choice.e_Null:
            obj = NullLoc.__new__(NullLoc)
        elif kind == CSeq_loc_choice.e_Empty:
            obj = EmptySeqLoc.__new__(EmptySeqLoc)
        elif kind == CSeq_loc_choice.e_Whole:
            obj = WholeSeqLoc.__new__(WholeSeqLoc)
        else:
            raise NotImplementedError(f"{kind!r}")

        obj._loc = ref
        return obj


cdef class NullLoc(SeqLoc):
    """A region of unknown length for which no data exists.
    """

cdef class EmptySeqLoc(SeqLoc):
    """A gap of unknown size inside an alignment.
    """

cdef class WholeSeqLoc(SeqLoc):
    """A reference to an entire `BioSeq`.
    """

    def __init__(self, SeqId sequence_id not None):
        """__init__(self, id)\n--\n

        Create a new location referencing the given sequence.

        Arguments:
            sequence_id (`~pyncbitk.objects.seqid.SeqId`): The identifier
                of the sequence being referenced. 

        """
        cdef CSeq_loc* loc = new CSeq_loc()
        self._loc.Reset(loc)
        loc.Select(CSeq_loc_choice.e_Whole)
        loc.SetWhole(sequence_id._ref.GetObject())

    def __reduce__(self):
        return type(self), (self.sequence_id,)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.sequence_id!r})"

    def __rich_repr__(self):
        yield self.sequence_id

    @property
    def sequence_id(self):
        """`~pyncbitk.objects.seqid.SeqId`: The identifier of the sequence.
        """
        id_ = CRef[CSeq_id](&self._loc.GetNonNullPointer().GetWholeMut())
        return SeqId._wrap(id_)


cdef class SeqIntervalLoc(SeqLoc):
    """A reference to an interval on a `BioSeq`.

    An interval is a single continuous region of defined length on a `BioSeq`.

    """
    
    # TODO: Handle strand.
    def __init__(self, SeqId sequence_id, TSeqPos start, TSeqPos stop, object strand = None):
        if start > stop:
            raise ValueError(f"interval limits in invalid order: {start} > {stop}")
        cdef CSeq_interval* inter = new CSeq_interval(sequence_id._ref.GetObject(), start, stop)
        cdef CSeq_loc* loc = new CSeq_loc()
        loc.Select(CSeq_loc_choice.e_Int)
        loc.SetInt(inter[0])
        self._loc.Reset(loc)

    # TODO: handle strand 
    def __reduce__(self):
        return type(self), (self.sequence_id, self.start, self.stop)

    # TODO: handle strand 
    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.sequence_id!r}, start={self.start!r}, stop={self.stop!r})"

    # TODO: handle strand 
    def __rich_repr__(self):
        yield self.sequence_id
        yield "start", self.start
        yield "stop", self.stop

    @property
    def sequence_id(self):
        """`~pyncbitk.objects.seqid.SeqId`: The identifier of the sequence.
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
        """`int`: The end of the sequence interval (inclusive).
        """
        # FIXME: use exclusive indexing like in Python?
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
