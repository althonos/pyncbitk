from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.string cimport string

from ...corelib.ncbitype cimport Uint4
from ...corelib.ncbimisc cimport TSeqPos
from ..seqcode.seq_code_type cimport ESeq_code_type
from .seq_data cimport CSeq_data, E_Choice as CSeq_data_choice

cdef extern from "objects/seq/seqport_util.hpp" namespace "ncbi::objects::CSeqportUtil" nogil:
    ctypedef unsigned int TIndex
    ctypedef pair[TIndex, TIndex] TPair

cdef extern from "objects/seq/seqport_util.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeqportUtil:

        @staticmethod
        TSeqPos Convert(
            const CSeq_data&       in_seq,
            CSeq_data*             out_seq,
            CSeq_data_choice       to_code,
            # TSeqPos                uBeginIdx = 0,
            # TSeqPos                uLength   = 0,
            # bool                   bAmbig    = false,
            # Uint4                  seed      = 17734276
        ) except +

        @staticmethod 
        TSeqPos ConvertWithBlastAmbig(
            const CSeq_data& in_seq,
            CSeq_data*       out_seq,
            TSeqPos          uBeginIdx,
            TSeqPos          uLength,
            TSeqPos          total_length,
            TSeqPos*         out_seq_length,
            vector[Uint4]*   blast_ambig
        ) except +
        
        @staticmethod
        TSeqPos Pack(
            CSeq_data*   in_seq,                        
            # TSeqPos uLength = ncbi::numeric_limits<TSeqPos>::max()
        ) noexcept

        @staticmethod
        bool FastValidate(
            const CSeq_data&   in_seq,
            # TSeqPos            uBeginIdx = 0,
            # TSeqPos            uLength   = 0
        ) noexcept

        @staticmethod
        void Validate(
            const CSeq_data&   in_seq,
            vector[TSeqPos]*   badIdx,
            # TSeqPos            uBeginIdx = 0,
            # TSeqPos            uLength   = 0
        ) noexcept

        @staticmethod
        TSeqPos GetAmbigs(
            const CSeq_data&    in_seq,
            CSeq_data*          out_seq,
            vector[TSeqPos]*    out_indices,
            CSeq_data_choice    to_code = CSeq_data_choice.e_Ncbi2na,
            # TSeqPos             uBeginIdx = 0,
            # TSeqPos             uLength   = 0
        ) noexcept

        @staticmethod
        TSeqPos GetCopy(
            const CSeq_data&   in_seq,
            CSeq_data*         out_seq,
            # TSeqPos            uBeginIdx = 0,
            # TSeqPos            uLength   = 0
        ) except +

        @staticmethod 
        TSeqPos Keep(
            CSeq_data*   in_seq,
            # TSeqPos      uBeginIdx = 0,
            # TSeqPos      uLength   = 0
        ) except +

        @staticmethod
        TSeqPos Append(
            CSeq_data*         out_seq,
            const CSeq_data&   in_seq1,
            TSeqPos            uBeginIdx1,
            TSeqPos            uLength1,
            const CSeq_data&   in_seq2,
            TSeqPos            uBeginIdx2,
            TSeqPos            uLength2
        ) except +

        @staticmethod 
        TSeqPos Complement(
            CSeq_data*   in_seq,
            # TSeqPos      uBeginIdx = 0,
            # TSeqPos      uLength   = 0
        ) except +

        @staticmethod 
        TSeqPos Complement(
            const CSeq_data&   in_seq,
            CSeq_data*         out_seq,
            # TSeqPos            uBeginIdx = 0,
            # TSeqPos            uLength   = 0
        ) except +

        @staticmethod
        TSeqPos Reverse(
            CSeq_data*   in_seq,
            # TSeqPos      uBeginIdx = 0,
            # TSeqPos      uLength   = 0
        ) except +

        @staticmethod 
        TSeqPos Reverse(
            const CSeq_data&   in_seq,
            CSeq_data*         out_seq,
            # TSeqPos            uBeginIdx = 0,
            # TSeqPos            uLength   = 0
        ) except +

        @staticmethod
        TSeqPos ReverseComplement(
            CSeq_data*   in_seq,
            # TSeqPos      uBeginIdx = 0,
            # TSeqPos      uLength   = 0
        ) except +

        @staticmethod 
        TSeqPos ReverseComplement(
            const CSeq_data&   in_seq,
            CSeq_data*         out_seq,
            # TSeqPos            uBeginIdx = 0,
            # TSeqPos            uLength   = 0
        ) except +
                                            
        @staticmethod
        const string& GetIupacaa3(TIndex ncbistdaa) except +
        
        @staticmethod
        bool IsCodeAvailable(CSeq_data_choice code_type) except +
        @staticmethod
        bool IsCodeAvailable(ESeq_code_type code_type) except +
        
        @staticmethod 
        TPair GetCodeIndexFromTo(CSeq_data_choice code_type) except +
        @staticmethod
        TPair GetCodeIndexFromTo(ESeq_code_type code_type) except +
            
        @staticmethod
        const string& GetCode(CSeq_data_choice code_type, TIndex idx) except +
        @staticmethod
        const string& GetCode(ESeq_code_type code_type, TIndex idx) except +
                            
        @staticmethod
        const string& GetName(CSeq_data_choice code_type, TIndex idx) except +
        @staticmethod
        const string& GetName(ESeq_code_type code_type, TIndex idx) except +

        @staticmethod
        TIndex GetIndex(CSeq_data_choice code_type, const string& code) except +
        @staticmethod
        TIndex GetIndex(ESeq_code_type code_type, const string& code) except +
        
        @staticmethod
        TIndex GetIndexComplement(CSeq_data_choice code_type, TIndex idx) except +
        @staticmethod
        TIndex GetIndexComplement(ESeq_code_type code_type, TIndex idx) except +
        
        @staticmethod
        TIndex GetMapToIndex(CSeq_data_choice from_type, CSeq_data_choice to_type, TIndex from_idx) except +
        @staticmethod
        TIndex GetMapToIndex(ESeq_code_type from_type, ESeq_code_type to_type, TIndex from_idx) except +