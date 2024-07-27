from libcpp cimport bool
from libcpp.set cimport set
from libcpp.string cimport string

from ....corelib.ncbiobj cimport CObject, CRef
from ....corelib.ncbitype cimport Int8
from .blast_types cimport EProgram
from .blast_options cimport EAPILocality, CBlastOptions


cdef extern from "algo/blast/api/blast_options_handle.hpp" namespace "ncbi::blast::CBlastOptionsFactory" nogil:

    enum ETaskSets:
        eNuclNucl
        eProtProt
        eMapping
        eAll


cdef extern from "algo/blast/api/blast_options_handle.hpp" namespace "ncbi::blast" nogil:

    cppclass CBlastOptionsFactory:
        @staticmethod
        CBlastOptionsHandle* Create(EProgram program)
        @staticmethod
        CBlastOptionsHandle* Create(EProgram program, EAPILocality locality)

        @staticmethod
        CBlastOptionsHandle* Create(string task)
        @staticmethod
        CBlastOptionsHandle* Create(string task, EAPILocality locality)

        @staticmethod
        set[string] GetTasks()
        @staticmethod
        set[string] GetTasks(ETaskSets choice)

    cppclass CBlastOptionsHandle(CObject):
        CBlastOptionsHandle(EAPILocality locality)
        bool Validate() except +

        const CBlastOptions& GetOptions() const
        CBlastOptions& GetOptionsRw "SetOptions" () noexcept

        void SetDefaults()
        void DoneDefaults()

        # Initial Word options
        int GetWindowSize() const
        void SetWindowSize(int ws)
        int GetOffDiagonalRange() const
        void SetOffDiagonalRange(int r)
        
        # Query setup options
        void ClearFilterOptions() 
        char* GetFilterString() const
        void SetFilterString(const char* f)
        void SetFilterString(const char* f, bool clear)

        bool GetMaskAtHash() const
        void SetMaskAtHash(bool m = true)

        double GetGapXDropoff() const
        void SetGapXDropoff(double x)

        double GetGapTrigger() const
        void SetGapTrigger(double g) 

        double GetGapXDropoffFinal() const
        void SetGapXDropoffFinal(double x)

        int GetHitlistSize() const
        void SetHitlistSize(int s)

        int GetMaxNumHspPerSequence() const
        void SetMaxNumHspPerSequence(int m)

        int GetMaxHspsPerSubject() const
        void SetMaxHspsPerSubject(int m)

        double GetEvalueThreshold() const 
        void SetEvalueThreshold(double eval_)
        int GetCutoffScore() const 
        void SetCutoffScore(int s)

        double GetPercentIdentity() const
        void SetPercentIdentity(double p)

        double GetQueryCovHspPerc() const
        void SetQueryCovHspPerc(double p)

        int GetMinDiagSeparation() const
        void SetMinDiagSeparation(int d)

        bool GetGappedMode() const
        void SetGappedMode()
        void SetGappedMode(bool m)

        int GetCullingLimit() const
        void SetCullingLimit(int s)

        int GetMaskLevel() const
        void SetMaskLevel(int ml)

        bool GetComplexityAdjMode() const
        void SetComplexityAdjMode()
        void SetComplexityAdjMode(bool m)

        double GetLowScorePerc() const
        void SetLowScorePerc(double p)

        Int8 GetDbLength() const 
        void SetDbLength(Int8 len)

        unsigned int GetDbSeqNum() const 
        void SetDbSeqNum(unsigned int num)

        Int8 GetEffectiveSearchSpace() const 
        void SetEffectiveSearchSpace(Int8 eff)
