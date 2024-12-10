from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..corelib.ncbiobj cimport CObject, CRef
from ..objects.seqloc.seq_id cimport CSeq_id
from ..objects.seqloc.seq_loc cimport CSeq_loc
from ..objects.seq.bioseq cimport CBioseq
from .object_manager cimport CObjectManager
from .bioseq_handle cimport CBioseq_Handle

cdef extern from "objmgr/scope.hpp" namespace "ncbi::objects::CScope" nogil:
    ctypedef int TPriority
    enum EPriority:
        kPriority_Default = -1
        kPriority_NotSet = -1
   
    enum EGetBioseqFlag:
        eGetBioseq_Resolved
        eGetBioseq_Loaded
        eGetBioseq_All

    # ctypedef vector[CSeq_id_Handle] TIds
    # ctypedef vector[CBioseq_Handle] TBioseqHandles
    # ctypedef vector[CTSE_Handle]    TCDD_Entries

    enum EMissing:
        eMissing_Throw
        eMissing_Null
        eMissing_Default

    # ctypedef CBlobIdKey TBlobId

    enum EExist:
        eExist_Throw
        eExist_Get
        eExist_Default

    enum EForceLoad:
        eNoForceLoad
        eForceLoad

    enum EGetFlags:
        fForceLoad
        fThrowOnMissingSequence
        fThrowOnMissingData
        fThrowOnMissing
        fDoNotRecalculate
    ctypedef int TGetFlags

    bool Exists(const CSeq_id& id)
    # bool Exists(const CSeq_id_Handle& id)

    # TIds GetIds(const CSeq_id&        id )
    # TIds GetIds(const CSeq_id&        id , TGetFlags flags)
    # TIds GetIds(const CSeq_id_Handle& idh)
    # TIds GetIds(const CSeq_id_Handle& idh, TGetFlags flags)


cdef extern from "objmgr/scope.hpp" namespace "ncbi::objects" nogil:
    
    cppclass CScope(CObject):
        CScope(CObjectManager& objmgr)

        CObjectManager& GetObjectManager() except +

        CBioseq_Handle GetBioseqHandle(const CSeq_id& id) except +
        # CBioseq_Handle GetBioseqHandle(const CSeq_id_Handle& id) except +
        CBioseq_Handle GetBioseqHandle(const CSeq_loc& loc) except +
        
        # CBioseq_Handle GetBioseqHandle(const CSeq_id& id, EGetBioseqFlag get_flag)
        # CBioseq_Handle GetBioseqHandle(const CSeq_id_Handle& id, EGetBioseqFlag get_flag)

        # bool IsSameBioseq(const CSeq_id_Handle& id1, const CSeq_id_Handle& id2, EGetBioseqFlag get_flag)

        # TBioseqHandles GetBioseqHandles(const TIds& ids)

        # TCDD_Entries GetCDDAnnots(const TIds& idhs)
        # TCDD_Entries GetCDDAnnots(const TBioseqHandles& bhs)

        # CSeq_entry_Handle GetSeq_entryHandle(
        #     CDataLoader* loader,
        #     const TBlobId& blob_id,
        #     EMissing action = eMissing_Default
        # )

        # CSeq_entry_Handle AddTopLevelSeqEntry(
        #     CSeq_entry& top_entry,
        # )
        # CSeq_entry_Handle AddTopLevelSeqEntry(
        #     CSeq_entry& top_entry,
        #     TPriority pri,
        # )
        # CSeq_entry_Handle AddTopLevelSeqEntry(
        #     CSeq_entry& top_entry,
        #     TPriority pri,
        #     EExist action
        # )
        # CSeq_entry_Handle AddTopLevelSeqEntry(
        #     const CSeq_entry& top_entry,
        # )
        # CSeq_entry_Handle AddTopLevelSeqEntry(
        #     const CSeq_entry& top_entry,
        #     TPriority pri,
        # )
        # CSeq_entry_Handle AddTopLevelSeqEntry(
        #     const CSeq_entry& top_entry,
        #     TPriority pri,
        #     EExist action,
        # )

        CBioseq_Handle AddBioseq(CBioseq& bioseq) except +
        # CBioseq_Handle AddBioseq(CBioseq& bioseq, TPriority pri)
        # CBioseq_Handle AddBioseq(CBioseq& bioseq, TPriority pri, EExist action)
        # CBioseq_Handle AddBioseq(const CBioseq& bioseq)
        # CBioseq_Handle AddBioseq(const CBioseq& bioseq, TPriority pri)
        # CBioseq_Handle AddBioseq(const CBioseq& bioseq, TPriority pri, EExist action)


        # Check existence of sequence with this id
        bool Exists(const CSeq_id&        id) except +
        # bool Exists(const CSeq_id_Handle& id)