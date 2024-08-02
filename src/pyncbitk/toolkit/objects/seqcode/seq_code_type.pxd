cdef extern from "objects/seqcode/Seq_code_type_.hpp" namespace "ncbi::objects" nogil:

    enum ESeq_code_type:
        eSeq_code_type_iupacna 
        eSeq_code_type_iupacaa 
        eSeq_code_type_ncbi2na  
        eSeq_code_type_ncbi4na  
        eSeq_code_type_ncbi8na  
        eSeq_code_type_ncbipna  
        eSeq_code_type_ncbi8aa  
        eSeq_code_type_ncbieaa  
        eSeq_code_type_ncbipaa  
        eSeq_code_type_iupacaa3 
        eSeq_code_type_ncbistdaa