# cython: language_level=3, linetrace=True, binding=True

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seqset.bioseq_set cimport CBioseq_set
from ..toolkit.objects.seqset.seq_entry cimport CSeq_entry

from ..serial cimport Serial

# --- BioSeqSet ----------------------------------------------------------------

cdef class BioSeqSet:
    cdef CRef[CBioseq_set] _ref

# --- Entry --------------------------------------------------------------------

cdef class Entry(Serial):
    cdef CRef[CSeq_entry] _ref

    @staticmethod
    cdef Entry _wrap(CRef[CSeq_entry] ref)

cdef class SeqEntry(Entry):
    pass

cdef class SetEntry(Entry):
    pass


