from libcpp cimport bool
from libcpp.list cimport list

from ...corelib.ncbiobj cimport CObject, CRef, CConstRef
from ...serial.serialbase cimport CSerialObject
from .seq_entry cimport CSeq_entry
from ..general.object_id cimport CObject_id


cdef extern from "objects/seqset/Bioseq_set_.hpp" namespace "ncbi::objects::CBioseq_set_Base" nogil:

    enum EClass:
        eClass_not_set          
        eClass_nuc_prot         
        eClass_segset           
        eClass_conset           
        eClass_parts            
        eClass_gibb             
        eClass_gi               
        eClass_genbank          
        eClass_pir              
        eClass_pub_set          
        eClass_equiv            
        eClass_swissprot        
        eClass_pdb_entry        
        eClass_mut_set          
        eClass_pop_set          
        eClass_phy_set          
        eClass_eco_set          
        eClass_gen_prod_set     
        eClass_wgs_set          
        eClass_named_annot      
        eClass_named_annot_prod 
        eClass_read_set         
        eClass_paired_end_reads 
        eClass_small_genome_set 

    ctypedef CObject_id TId
    # ctypedef CDbtag TColl
    # ctypedef int TLevel
    # ctypedef EClass TClass
    # ctypedef string TRelease
    # ctypedef CDate TDate
    # ctypedef CSeq_descr TDescr
    ctypedef list[CRef[CSeq_entry]] TSeq_set
    # ctypedef list[CRef[CSeq_annot]] TAnnot
    

cdef extern from "objects/seqset/Bioseq_set_.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_set_Base(CSerialObject):
        CBioseq_set_Base()

        bool IsSetId() const
        bool CanGetId() const
        void ResetId()
        const TId& GetId() const
        void SetId(TId& value)
        TId& GetIdMut "SetId"()

        bool IsSetSeq_set() const
        bool CanGetSeq_set() const
        void ResetSeq_set()
        const TSeq_set& GetSeq_set() const
        TSeq_set& GetSeq_setMut "SetSeq_set"()


cdef extern from "objects/seqset/Bioseq_set.hpp" namespace "ncbi::objects" nogil:

    cppclass CBioseq_set(CBioseq_set_Base):
        CBioseq_set()