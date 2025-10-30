# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from .toolkit.objmgr.object_manager cimport CObjectManager
from .toolkit.objmgr.scope cimport CScope
from .toolkit.corelib.ncbiobj cimport CObject, CConstRef, CRef
from .toolkit.objects.seqloc.seq_id cimport CSeq_id
from .toolkit.objects.seq.seq_inst cimport CSeq_inst
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objects.seq.seq_descr cimport CSeq_descr

from .objects.seq cimport BioSeq
from .objects.seqinst cimport SeqInst
from .objects.seqloc cimport SeqLoc
from .objects.seqid cimport SeqId
from .objects.seqdesc cimport SeqDescSet


cdef class BioSeqHandle:

    def __call__(self):
        """Resolve the handle to load the complete `BioSeq`.
        """
        cdef CConstRef[CBioseq] ccref = self._handle.GetCompleteBioseq()
        cdef BioSeq bioseq = BioSeq.__new__(BioSeq)
        bioseq._ref.Reset(<CBioseq*> ccref.GetNonNullPointer())
        return bioseq

    @property
    def scope(self):
        cdef Scope scope = Scope.__new__(Scope)
        scope._scope.Reset(&self._handle.GetScope())
        return scope

    @property
    def id(self):
        """`SeqId` or `None`: The identifier for this handle, if any.
        """
        cdef CConstRef[CSeq_id] id_ = self._handle.GetInitialSeqIdOrNull()
        if id_.IsNull():
            return None
        cdef CRef[CSeq_id] ref = CRef[CSeq_id](&id_.GetObject())
        return SeqId._wrap(ref)

    @property
    def ids(self):
        """`list` of `~pyncbitk.objects.seqloc.SeqId`: The sequence identifiers.
        """
        raise NotImplementedError("BioSeqHandle.ids")

    @property
    def instance(self):
        """`~pyncbitk.objects.seqinst.SeqInst`: The sequence instance.
        """
        assert self._handle.CanGetInst()
        cdef CSeq_inst* inst = <CSeq_inst*> &self._handle.GetInst()
        return SeqInst._wrap(CRef[CSeq_inst](inst))

    @property
    def length(self):
        """`int`: The length of the sequence instance.
        """
        return self._handle.GetInst_Length()

    @property
    def descriptions(self):
        """`~pyncbitk.objects.seqdesc.SeqDescSet`: The sequence descriptions.
        """
        assert self._handle.CanGetDescr()
        cdef CSeq_descr* desc = <CSeq_descr*> &self._handle.GetDescr()
        return SeqDescSet._wrap(CRef[CSeq_descr](desc))



cdef class ObjectManager:
    """The global object manager.
    """

    def __init__(self):
        self._mgr = CObjectManager.GetInstance()

    cpdef Scope scope(self):
        """Create a new scope.
        """
        return Scope(self)

    def register_data_loader(self, str name not None):
        cdef bytes _name = name.encode()
        self._mgr.GetObject().RegisterDataLoader(NULL, <string> _name)

    def get_registered_names(self):
        cdef vector[string] names
        self._mgr.GetObject().GetRegisteredNames(names)
        return list(names)


cdef class Scope:
    """A handler for managing objects within a given scope.

    This class works as a context manager to support adding data to the 
    object manager for a given scope.

    """

    def __init__(self, ObjectManager manager):
        self._scope.Reset(new CScope(manager._mgr.GetObject()))

    def __enter__(self):
        return self

    def __exit__(self, exc_ty, exc_val, traceback):
        self.close()
        return False

    def __contains__(self, object key):
        if not isinstance(key, SeqId):
            return False
        cdef SeqId id_ = key
        return self._scope.GetObject().Exists(id_._ref.GetObject())

    def __getitem__(self, SeqId key):
        cdef BioSeqHandle handle = BioSeqHandle.__new__(BioSeqHandle)
        handle._handle = self._scope.GetObject().GetBioseqHandle(key._ref.GetObject())
        return handle

    def close(self):
        """Close the scope and release the associated data.
        """
        self._scope.ReleaseOrNull()

    cpdef BioSeqHandle add_bioseq(self, BioSeq seq) except *:
        cdef CBioseq_Handle handle = self._scope.GetObject().AddBioseq(seq._ref.GetObject())
        cdef BioSeqHandle obj = BioSeqHandle.__new__(BioSeqHandle) 
        obj._handle = handle
        return obj
