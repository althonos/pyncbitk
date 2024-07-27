from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector


from ....corelib.ncbiobj cimport CObject, CRef
from ....corelib.ncbitype cimport Uint4 
from ....objects.seqalign.seq_align_set cimport CSeq_align_set

cdef extern from "algo/blast/api/sseqloc.hpp" namespace "ncbi::blast" nogil:

    enum EProgram:
        eBlastNotSet
        eBlastn
        eBlastp
        eBlastx
        eTblastn
        eTblastx
        eRPSBlast
        eRPSTblastn
        eMegablast
        eDiscMegablast
        ePSIBlast
        ePSITblastn
        ePHIBlastp
        ePHIBlastn
        eDeltaBlast
        eVecScreen
        eMapper
        eKBlastp
        eBlastProgramMax

    string EProgramToTaskName(EProgram p) except +
    EProgram ProgramNameToEnum(const string& program_name) except +
    void ThrowIfInvalidTask(const string& task) except +
    # EBlastProgramType EProgramToEBlastProgramType(EProgram p)

    cppclass CSearchMessage(CObject):
        # CSearchMessage(EBlastSeverity severity, int error_id, const string& message)
        CSearchMessage()

        # EBlastSeverity GetSeverity() const
        # void SetSeverity(EBlastSeverity sev)
        # string GetSeverityString() const
        # @staticmethod
        # string GetSeverityString(EBlastSeverity severity) const

        int GetErrorId() const
        string& SetMessage()
        string GetMessage() const
        string GetMessage(bool withSeverity) const

    cppclass TQueryMessages(vector[CRef[CSearchMessage]]):
        void SetQueryId(const string& id)
        string GetQueryId() const
        void Combine(const TQueryMessages& other)

    cppclass TSearchMessages(vector[CRef[TQueryMessages]]):
        # void AddMessageAllQueries(EBlastSeverity severity, int error_id, const string& message)
        bool HasMessages() const
        string ToString() const
        void Combine(const TSearchMessages& other_msgs)
        void RemoveDuplicates()

    enum EResultType:
        eDatabaseSearch
        eSequenceComparison

    ctypedef vector[CRef[CSeq_align_set]] TSeqAlignVector