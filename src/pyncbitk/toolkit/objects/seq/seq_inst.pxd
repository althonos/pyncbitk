from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.ncbimisc cimport TSeqPos
from ...serial.serialbase cimport CSerialObject
from ...serial.serialdef cimport ESerialRecursionMode
from .seq_data cimport CSeq_data
from .seq_ext cimport CSeq_ext

cdef extern from "objects/seq/Seq_inst_.hpp" namespace "ncbi::objects::CSeq_inst_Base" nogil:
    
    enum ERepr:
        eRepr_not_set
        eRepr_virtual
        eRepr_raw    
        eRepr_seg    
        eRepr_const 
        eRepr_ref   
        eRepr_consen
        eRepr_map   
        eRepr_delta 
        eRepr_other 

    enum EMol:
        eMol_not_set
        eMol_dna    
        eMol_rna    
        eMol_aa     
        eMol_na     
        eMol_other  

    enum ETopology:
        eTopology_not_set 
        eTopology_linear  
        eTopology_circular 
        eTopology_tandem   
        eTopology_other   

    enum EStrand:
        eStrand_not_set
        eStrand_ss     
        eStrand_ds      
        eStrand_mixed 
        eStrand_other   

    enum E_memberIndex:
        e__allMandatory
        e_repr
        e_mol
        e_length
        e_fuzz
        e_topology
        e_strand
        e_seq_data
        e_ext
        e_hist

    ctypedef ERepr TRepr
    ctypedef EMol TMol
    ctypedef TSeqPos TLength
    # ctypedef CInt_fuzz TFuzz
    ctypedef ETopology TTopology
    ctypedef EStrand TStrand
    ctypedef CSeq_data TSeq_data
    ctypedef CSeq_ext TExt
    # ctypedef CSeq_hist THist

cdef extern from "objects/seq/Seq_inst_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_inst_Base(CSerialObject):
        CSeq_inst_Base()

        # mandatory
        bool IsSetRepr() const
        bool CanGetRepr() const
        void ResetRepr()
        TRepr GetRepr() const
        void SetRepr(TRepr value)
        TRepr& GetReprMut "SetRepr"()

        # mandatory
        bool IsSetMol() const
        bool CanGetMol() const
        void ResetMol()
        TMol GetMol() const
        void SetMol(TMol value)
        TMol& GetMolMut "SetMol"()

        # optional
        bool IsSetLength() const
        bool CanGetLength() const
        void ResetLength()
        TLength GetLength() const
        void SetLength(TLength value)
        TLength& GetLengthMut "SetLength"()

        # # optional
        # bool IsSetFuzz() const
        # bool CanGetFuzz() const
        # void ResetFuzz()
        # const TFuzz& GetFuzz() const
        # void SetFuzz(TFuzz& value)
        # TFuzz& GetFuzzRw "SetFuzz"()

        # optional with default eTopology_linear
        bool IsSetTopology() const
        bool CanGetTopology() const
        void ResetTopology()
        void SetDefaultTopology()
        TTopology GetTopology() const
        void SetTopology(TTopology value)
        TTopology& GetTopologyMut "SetTopology"()

        # optional
        bool IsSetStrand() const
        bool CanGetStrand() const
        void ResetStrand()
        TStrand GetStrand() const
        void SetStrand(TStrand value)
        TStrand& GetStrandMut "SetStrand" ()

        # optional
        bool IsSetSeq_data() const
        bool CanGetSeq_data() const
        void ResetSeq_data()
        const TSeq_data& GetSeq_data() const
        void SetSeq_data(TSeq_data& value)
        TSeq_data& GetSeq_dataMut "SetSeq_data" ()

        # optional
        bool IsSetExt() const
        bool CanGetExt() const
        void ResetExt()
        const TExt& GetExt() const
        void SetExt(TExt& value)
        TExt& GetExtMut "SetExt"()

        # # optional
        # bool IsSetHist() const
        # bool CanGetHist() const
        # void ResetHist()
        # const THist& GetHist() const
        # void SetHist(THist& value)
        # THist& SetHist()

        void Reset()

cdef extern from "objects/seq/Seq_inst.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_inst(CSeq_inst_Base):
        CSeq_inst()

        bool IsNa()
        bool IsAa()

        @staticmethod
        string GetMoleculeClass(EMol mol)

        bool ConvertDeltaToRaw()


