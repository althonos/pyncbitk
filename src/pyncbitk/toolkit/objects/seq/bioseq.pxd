from libcpp cimport bool
from libcpp.string cimport string
from libcpp.list cimport list as cpplist

from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from ..seqloc.seq_id cimport CSeq_id
from .seq_inst cimport CSeq_inst
from .seq_descr cimport CSeq_descr

cdef extern from "objects/seq/Bioseq_.hpp" namespace "ncbi::objects::CBioseq_Base" nogil:

    ctypedef cpplist[CRef[CSeq_id]] TId
    ctypedef CSeq_descr TDescr
    ctypedef CSeq_inst TInst
    # ctypedef cpplist[CRef[CSeq_annot]] TAnnot


cdef extern from "objects/seq/Bioseq_.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_Base(CSerialObject):
        CBioseq_Base()

        # mandatory
        bool IsSetId() const
        bool CanGetId() const
        void ResetId()
        const TId& GetId() const
        TId& SetId()

        # mandatory
        bool IsSetDescr() const
        bool CanGetDescr() const
        void ResetDescr()
        const TDescr& GetDescr() except +
        void SetDescr(TDescr& value) except +
        TDescr& GetDescrMut "SetDescr" () except +

        # mandatory
        bool IsSetInst() const
        bool CanGetInst() const
        void ResetInst()
        const TInst& GetInst() except +
        void SetInst(TInst& value) except +
        TInst& GetInstMut "SetInst" () except +

        # optional
        # bool IsSetAnnot() const
        # bool CanGetAnnot() const
        # void ResetAnnot()
        # const TAnnot& GetAnnot() const
        # TAnnot& SetAnnot()

        void Reset()


cdef extern from "objects/seq/Bioseq.hpp" namespace "ncbi::objects::CBioseq" nogil:

    enum ELabelType:
        eType
        eContent
        eBoth

cdef extern from "objects/seq/Bioseq.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq(CBioseq_Base):
        CBioseq()

        void GetLabel(string* label, ELabelType type, bool worst = false) const

        const CSeq_id* GetFirstId() const
        const CSeq_id* GetNonLocalId() const
        const CSeq_id* GetLocalId() const

        bool IsNa() const
        bool IsAa() const
        
