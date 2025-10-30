# cython: language_level=3

from libcpp.list cimport list as cpplist

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seqset.seq_entry cimport CSeq_entry, E_Choice as CSeq_entry_choice
from ..toolkit.serial.serialbase cimport CSerialObject, EResetVariant
from ..toolkit.objects.seq.bioseq cimport CBioseq
from ..toolkit.objects.seqset.bioseq_set cimport CBioseq_set

from .general cimport ObjectId
from .seqset cimport Entry
from .seq cimport BioSeq

# --- BioSeqSet ----------------------------------------------------------------

cdef class BioSeqSet(Serial):
    """A set of biological sequence.

    A `BioSeqSet` is a set that stores either other sequences (as `BioSeq`
    objects), or recursively other sequence sets (as `BioSeqSet`),
    allowing to create tree of sequences.

    """

    @staticmethod
    cdef BioSeqSet _wrap(CRef[CBioseq_set] ref):
        cdef BioSeqSet seqset = BioSeqSet.__new__(BioSeqSet)
        seqset._ref = ref
        return seqset

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    # TODO: subclasses

    def __init__(
        self,
        object items=(),
        *,
        ObjectId id = None
    ):
        cdef Entry        entry
        cdef CBioseq_set* seqset = new CBioseq_set()

        try:
            if id is not None:
                seqset.SetId(id._ref.GetObject())
            for item in items:
                if isinstance(item, Entry):
                    entry = item
                elif isinstance(item, BioSeq):
                    entry = SeqEntry(item)
                elif isinstance(item, BioSeqSet):
                    entry = SetEntry(item)
                else:
                    ty = item.__class__.__name__
                    raise TypeError(f"expected Entry, BioSeq or BioSeqSet, got {ty}")
                seqset.GetSeq_setMut().push_back(entry._ref)
        finally:
            self._ref.Reset(seqset)

    def __len__(self):
        assert self._ref.GetObject().IsSetSeq_set()
        return self._ref.GetObject().GetSeq_set().size()

    def __iter__(self):
        cdef CRef[CSeq_entry]          item
        cdef cpplist[CRef[CSeq_entry]] items = self._ref.GetObject().GetSeq_setMut()
        for item in items:
            yield Entry._wrap(item)

    def __rich_repr__(self):
        yield list(self)

# --- Entry --------------------------------------------------------------------

cdef class Entry(Serial):
    """A sequence entry.
    """

    @staticmethod
    cdef Entry _wrap(CRef[CSeq_entry] ref):

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

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

cdef class SeqEntry(Entry):

    def __init__(self, BioSeq sequence):
        cdef CSeq_entry* entry = new CSeq_entry()
        entry.Select(CSeq_entry_choice.e_Seq, EResetVariant.eDoResetVariant)
        entry.SetSeq(sequence._ref.GetObject())
        self._ref.Reset(entry)

    @property
    def sequence(self):
        cdef CBioseq* bioseq = &self._ref.GetNonNullPointer().GetSeqMut()
        return BioSeq._wrap(CRef[CBioseq](bioseq))


cdef class SetEntry(Entry):

    def __init__(self, BioSeqSet set):
        cdef CSeq_entry* entry = new CSeq_entry()
        entry.Select(CSeq_entry_choice.e_Set, EResetVariant.eDoResetVariant)
        entry.SetSet(set._ref.GetObject())
        self._ref.Reset(entry)

    @property
    def set(self):
        cdef CBioseq_set* bioseq_set = &self._ref.GetNonNullPointer().GetSetMut()
        return BioSeqSet._wrap(CRef[CBioseq_set](bioseq_set))


