from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..corelib.ncbiobj cimport CObject, CRef, CConstRef
from ..objects.seqloc.seq_id cimport CSeq_id
from ..objects.seq.bioseq cimport CBioseq
from ..objects.seq.seq_inst cimport CSeq_inst, TLength as TInstLength
from ..objects.seq.seq_descr cimport CSeq_descr
# from .object_manager cimport CObjectManager
from .scope cimport CScope


cdef extern from "objmgr/bioseq_handle.hpp" namespace "ncbi::objects::CBioseq_Handle" nogil:

    ctypedef CBioseq TObject
    ctypedef CSeq_inst TInst
    ctypedef TInstLength TInst_Length
    ctypedef CSeq_descr TDescr


cdef extern from "objmgr/bioseq_handle.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_Handle:
        CBioseq_Handle()
        
        void Reset()
        CScope& GetScope() const

        CConstRef[CSeq_id] GetSeqId() const
        CConstRef[CSeq_id] GetInitialSeqIdOrNull() const

        CConstRef[CBioseq] GetCompleteBioseq() except +
        
        CConstRef[TObject] GetCompleteObject() except +
        CConstRef[TObject] GetObjectCore() except +

        # inst
        # typedef vector<CSeq_id_Handle> TId;
        # bool IsSetId(void) const;
        # bool CanGetId(void) const;
        # const TId& GetId(void) const;
        
        # descr
        bool IsSetDescr() const
        bool CanGetDescr() const
        const TDescr& GetDescr() except +
        
        # inst
        bool IsSetInst() const
        bool CanGetInst() const
        const TInst& GetInst() except +
        
        # // inst.repr
        # typedef TInst::TRepr TInst_Repr;
        # bool IsSetInst_Repr(void) const;
        # bool CanGetInst_Repr(void) const;
        # TInst_Repr GetInst_Repr(void) const;

        # // inst.mol
        # typedef TInst::TMol TInst_Mol;
        # bool IsSetInst_Mol(void) const;
        # bool CanGetInst_Mol(void) const;
        # TInst_Mol GetInst_Mol(void) const;

        # // inst.length
        bool IsSetInst_Length() const
        bool CanGetInst_Length() const
        TInst_Length GetInst_Length() except +

        # TSeqPos GetBioseqLength(void) const; // try to calculate it if not set
        # // inst.fuzz
        # typedef TInst::TFuzz TInst_Fuzz;
        # bool IsSetInst_Fuzz(void) const;
        # bool CanGetInst_Fuzz(void) const;
        # const TInst_Fuzz& GetInst_Fuzz(void) const;

        # // inst.topology
        # typedef TInst::TTopology TInst_Topology;
        # bool IsSetInst_Topology(void) const;
        # bool CanGetInst_Topology(void) const;
        # TInst_Topology GetInst_Topology(void) const;
        # // inst.strand
        # typedef TInst::TStrand TInst_Strand;
        # bool IsSetInst_Strand(void) const;
        # bool CanGetInst_Strand(void) const;
        # TInst_Strand GetInst_Strand(void) const;
        # // inst.seq-data
        # typedef TInst::TSeq_data TInst_Seq_data;
        # bool IsSetInst_Seq_data(void) const;
        # bool CanGetInst_Seq_data(void) const;
        # const TInst_Seq_data& GetInst_Seq_data(void) const;
        # // inst.ext
        # typedef TInst::TExt TInst_Ext;
        # bool IsSetInst_Ext(void) const;
        # bool CanGetInst_Ext(void) const;
        # const TInst_Ext& GetInst_Ext(void) const;
        # // inst.hist
        # typedef TInst::THist TInst_Hist;
        # bool IsSetInst_Hist(void) const;
        # bool CanGetInst_Hist(void) const;
        # const TInst_Hist& GetInst_Hist(void) const;
        # // annot
        # bool HasAnnots(void) const;

        # // Check sequence type
        # typedef CSeq_inst::TMol TMol;
        # TMol GetSequenceType(void) const;
        # bool IsProtein(void) const;
        # bool IsNucleotide(void) const;

        bool operator== (const CBioseq_Handle& h) const
        bool operator!= (const CBioseq_Handle& h) const
        bool operator<  (const CBioseq_Handle& h) const

        bool IsRemoved() const