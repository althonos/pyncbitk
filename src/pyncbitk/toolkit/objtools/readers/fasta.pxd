from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.ncbiobj cimport CRef
from ...objects.seqset.seq_entry cimport CSeq_entry
from .reader_base cimport CReaderBase

cdef extern from "objtools/readers/fasta.hpp" namespace "ncbi::objects::CFastaReader" nogil:

    ctypedef long TFlags
    
    enum EFlags:
        fAssumeNuc    
        fAssumeProt    
        fForceType      
        fNoParseID    
        fParseGaps    
        fOneSeq     
        fAllSeqIds 
        fNoSeqData    
        fRequireID    
        fDLOptional   
        fParseRawID   
        fSkipCheck    
        fNoSplit      
        fValidate     
        fUniqueIDs    
        fStrictGuess  
        fLaxGuess     
        fAddMods          
        fLetterGaps       
        fNoUserObjs   
        fBadModThrow  
        fUnknModThrow 
        fLeaveAsText  
        fQuickIDCheck 
        fUseIupacaa   
        fHyphensIgnoreAndWarn
        fDisableNoResidues   
        fDisableParseRange   
        fIgnoreMods          


cdef extern from "objtools/readers/fasta.hpp" namespace "ncbi::objects" nogil:

    cppclass CFastaReader(CReaderBase):
        # CFastaReader(ILineReader& reader, TFlags flags = 0, FIdCheck f_idcheck = CSeqIdCheck());
        # CFastaReader(CNcbiIstream& in,    TFlags flags = 0, FIdCheck f_idcheck = CSeqIdCheck());
        # CFastaReader(const string& path,  TFlags flags = 0, FIdCheck f_idcheck = CSeqIdCheck());
        CFastaReader(const string& path)
        CFastaReader(const string& path, TFlags flags)
        # CFastaReader(const string& path, TFlags flags)
        # CFastaReader(CReaderBase::TReaderFlags fBaseFlags,
        #             TFlags flags = 0,
        #             FIdCheck f_idcheck = CSeqIdCheck()
        #             );

        bool AtEOF() const
        CRef[CSeq_entry] ReadOneSeq() except +
        # virtual CRef<CSeq_entry> ReadOneSeq(ILineErrorListener* pMessageListener = nullptr);

