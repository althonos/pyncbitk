# cython: language_level=3, linetrace=True, binding=True

from libcpp.list cimport list as cpplist

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seqset.seq_entry cimport CSeq_entry, E_Choice as CSeq_entry_choice
from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.seq.bioseq cimport CBioseq

from .seqset cimport Entry
from .seq cimport BioSeq

# --- BioSeqSet ----------------------------------------------------------------

cdef class BioSeqSet:
    """A set of sequence entries.
    """

    # TODO: subclasses

    def __init__(self):
        pass

    def __len__(self):
        assert self._ref.GetObject().IsSetSeq_set()
        return self._ref.GetObject().GetSeq_set().size()

    def __iter__(self):
        cdef CRef[CSeq_entry]          item
        cdef cpplist[CRef[CSeq_entry]] items = self._ref.GetObject().GetSeq_setMut()
        for item in items:
            yield Entry._wrap(item)

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

    @property
    def seq(self):
        cdef CBioseq* bioseq = &self._ref.GetNonNullPointer().SetSeq()
        return BioSeq._wrap(CRef[CBioseq](bioseq))


cdef class SetEntry(Entry):
    pass


