# cython: language_level=3

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seqset.bioseq_set cimport CBioseq_set
from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.objects.seqset.seq_entry cimport CSeq_entry
from ..toolkit.objects.seqset.bioseq_set cimport CBioseq_set

from ..serial cimport Serial

# --- BioSeqSet ----------------------------------------------------------------

cdef class BioSeqSet(Serial):
    cdef CRef[CBioseq_set] _ref

    @staticmethod
    cdef BioSeqSet _wrap(CRef[CBioseq_set] ref)

    cdef CSerialObject* _serial(self)


# --- Entry --------------------------------------------------------------------

cdef class Entry(Serial):
    cdef CRef[CSeq_entry] _ref

    @staticmethod
    cdef Entry _wrap(CRef[CSeq_entry] ref)

cdef class SeqEntry(Entry):
    pass

cdef class SetEntry(Entry):
    pass


