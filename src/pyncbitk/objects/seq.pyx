# cython: language_level=3, linetrace=True, binding=True

from libcpp.list cimport list as cpplist

from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.objects.seq.bioseq cimport CBioseq
from ..toolkit.objects.seq.seq_inst cimport CSeq_inst
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id
from ..toolkit.serial.serialbase cimport CSerialObject

from ..serial cimport Serial
from .seqinst cimport SeqInst
from .seqloc cimport SeqId, LocalId

# --- BioSeq -------------------------------------------------------------------

cdef class BioSeq(Serial):
    """A biological sequence.
    """

    @staticmethod
    cdef BioSeq _wrap(CRef[CBioseq] ref):
        cdef BioSeq obj = BioSeq.__new__(BioSeq)
        obj._ref = ref
        return obj

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()

    def __init__(
        self,
        object ids not None,
        SeqInst instance not None,
    ):
        cdef SeqId    id_
        cdef CBioseq* obj = new CBioseq()

        obj.SetInst(instance._ref.GetObject())
        for id_ in ids:
            obj.SetId().push_back(id_._ref)
        if obj.GetId().size() == 0:
            raise ValueError("BioSeq must have at least one identifier")

        self._ref.Reset(obj)

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.ids!r}, {self.instance!r})"

    @property
    def ids(self):
        """`list` of `SeqId`: The identifiers of the sequence.
        """
        assert self._ref.GetNonNullPointer().IsSetId()  # mandatory

        cdef SeqId                  seqid
        cdef list                   ids   = []
        cdef cpplist[CRef[CSeq_id]] _ids  = self._ref.GetObject().GetId()

        for ref in _ids:
            ids.append(SeqId._wrap(ref))

        return ids

    @property
    def instance(self):
        """`SeqInst`: The actual sequence instance.
        """
        assert self._ref.GetNonNullPointer().IsSetInst()  # mandatory
        return SeqInst._wrap(CRef[CSeq_inst](&self._ref.GetObject().GetInstMut()))
