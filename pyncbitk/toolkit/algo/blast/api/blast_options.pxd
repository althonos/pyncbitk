from libcpp cimport bool

from ....corelib.ncbiobj cimport CObject, CRef
from .blast_types cimport EProgram


cdef extern from "algo/blast/api/blast_options.hpp" namespace "ncbi::blast::CBlastOptions" nogil:

    enum EAPILocality:
        eLocal
        eRemote
        eBoth

cdef extern from "algo/blast/api/blast_options.hpp" namespace "ncbi::blast" nogil:

    cppclass CBlastOptions(CObject):
        CBlastOptions()
        CBlastOptions(EAPILocality locality)

        CRef[CBlastOptions] Clone() const
        EAPILocality GetLocality() const

        bool Validate() const

        EProgram GetProgram() const
        void SetProgram(EProgram p)

        # EBlastProgramType GetProgramType() const
        bool IsIterativeSearch() const

        double GetWordThreshold() const
        void SetWordThreshold(double w)

        # ELookupTableType GetLookupTableType() const;
        # void SetLookupTableType(ELookupTableType type);

        # int GetWordSize() const
        # void SetWordSize(int ws)

        # Uint4 GetLookupTableStride() const
        # void SetLookupTableStride(Uint4 val)

        # bool GetLookupDbFilter(void) const
        # void SetLookupDbFilter(bool val)

        # Uint1 GetMaxDbWordCount(void) const
        # void SetMaxDbWordCount(Uint1 num)

        # unsigned char GetMBTemplateLength() const
        # void SetMBTemplateLength(unsigned char len)

        # unsigned char GetMBTemplateType() const
        # void SetMBTemplateType(unsigned char type)

        # void ClearFilterOptions()

        # NCBI_DEPRECATED char* GetFilterString() const
        # NCBI_DEPRECATED void SetFilterString(const char* f, bool clear = true)

        # bool GetMaskAtHash() const
        # void SetMaskAtHash(bool val = true)

        # bool GetDustFiltering() const
        # void SetDustFiltering(bool val = true)

        # int GetDustFilteringLevel() const
        # void SetDustFilteringLevel(int m)

        # int GetDustFilteringWindow() const
        # void SetDustFilteringWindow(int m)

        # int GetDustFilteringLinker() const
        # void SetDustFilteringLinker(int m)

        # bool GetSegFiltering() const
        # void SetSegFiltering(bool val = true)

        # int GetSegFilteringWindow() const
        # void SetSegFilteringWindow(int m)

        # double GetSegFilteringLocut() const
        # void SetSegFilteringLocut(double m)

        # double GetSegFilteringHicut() const
        # void SetSegFilteringHicut(double m)

        # bool GetRepeatFiltering() const;
        # void SetRepeatFiltering(bool val = true);

        # const char* GetRepeatFilteringDB() const;
        # void SetRepeatFilteringDB(const char* db);

        # int GetWindowMaskerTaxId() const;
        # void SetWindowMaskerTaxId(int taxid);

        # const char* GetWindowMaskerDatabase() const
        # void SetWindowMaskerDatabase(const char* db)

        # bool GetReadQualityFiltering() const
        # void SetReadQualityFiltering()
        # void SetReadQualityFiltering(bool val)

        # double GetReadMaxFractionAmbiguous() const
        # void SetReadMaxFractionAmbiguous(double val)

        # int GetReadMinDimerEntropy() const
        # void SetReadMinDimerEntropy(int val)


        # objects::ENa_strand GetStrandOption() const;
        # void SetStrandOption(objects::ENa_strand s);

        # int GetQueryGeneticCode() const;
        # void SetQueryGeneticCode(int gc);

        # int GetWindowSize() const;
        # void SetWindowSize(int w);

        # int GetOffDiagonalRange() const;
        # void SetOffDiagonalRange(int r);

        # double GetXDropoff() const;
        # void SetXDropoff(double x);

        # double GetGapXDropoff() const;
        # void SetGapXDropoff(double x);

        # double GetGapXDropoffFinal() const;
        # void SetGapXDropoffFinal(double x);

        # double GetGapTrigger() const;
        # void SetGapTrigger(double g);

        # EBlastPrelimGapExt GetGapExtnAlgorithm() const;
        # void SetGapExtnAlgorithm(EBlastPrelimGapExt a);

        # EBlastTbackExt GetGapTracebackAlgorithm() const;
        # void SetGapTracebackAlgorithm(EBlastTbackExt a);

        # ECompoAdjustModes GetCompositionBasedStats() const;
        # void SetCompositionBasedStats(ECompoAdjustModes mode);

        # bool GetSmithWatermanMode() const;
        # void SetSmithWatermanMode(bool m = true);

        # int GetUnifiedP() const;
        # void SetUnifiedP(int u = 0);

        # int GetMaxMismatches() const;
        # void SetMaxMismatches(int m);

        # int GetMismatchWindow() const;
        # void SetMismatchWindow(int w);

        # int GetHitlistSize() const;
        # void SetHitlistSize(int s);

        # int GetMaxNumHspPerSequence() const;
        # void SetMaxNumHspPerSequence(int m);

        # int GetMaxHspsPerSubject() const;
        # void SetMaxHspsPerSubject(int m);

        # int GetCullingLimit() const;
        # void SetCullingLimit(int s);

        # bool GetSubjectBestHit() const;
        # void SetSubjectBestHit();

        # double GetBestHitOverhang() const;
        # void SetBestHitOverhang(double overhang);
        # double GetBestHitScoreEdge() const;
        # void SetBestHitScoreEdge(double score_edge);

        # double GetEvalueThreshold() const;
        # void SetEvalueThreshold(double eval);

        # int GetCutoffScore() const;
        # void SetCutoffScore(int s);

        # vector<double> GetCutoffScoreCoeffs() const;
        # void SetCutoffScoreCoeffs(const vector<double>& c);

        # double GetPercentIdentity() const;
        # void SetPercentIdentity(double p);

        # int GetMaxEditDistance() const;
        # void SetMaxEditDistance(int e);

        # double GetQueryCovHspPerc() const;
        # void SetQueryCovHspPerc(double p);

        # int GetMinDiagSeparation() const;
        # void SetMinDiagSeparation(int d);

        # bool GetSumStatisticsMode() const;
        # void SetSumStatisticsMode(bool m = true);

        # int GetLongestIntronLength() const;
        # void SetLongestIntronLength(int l);

        # bool GetGappedMode() const;
        # void SetGappedMode(bool m = true);

        # int GetMaskLevel() const;
        # void SetMaskLevel(int s);

        # bool GetComplexityAdjMode() const;
        # void SetComplexityAdjMode(bool m = true);

        # double GetLowScorePerc() const;
        # void SetLowScorePerc(double p = 0.0);

        # bool GetPaired() const;
        # void SetPaired(bool p);

        # bool GetSpliceAlignments() const;
        # void SetSpliceAlignments(bool s);

        # const char* GetMatrixName() const;
        # void SetMatrixName(const char* matrix);

        # int GetMatchReward() const;
        # void SetMatchReward(int r);

        # int GetMismatchPenalty() const;
        # void SetMismatchPenalty(int p);

        # int GetGapOpeningCost() const;
        # void SetGapOpeningCost(int g);

        # int GetGapExtensionCost() const;
        # void SetGapExtensionCost(int e);

        # int GetFrameShiftPenalty() const;
        # void SetFrameShiftPenalty(int p);

        # bool GetOutOfFrameMode() const;
        # void SetOutOfFrameMode(bool m = true);

        # Int8 GetDbLength() const;
        # void SetDbLength(Int8 l);

        # unsigned int GetDbSeqNum() const;
        # void SetDbSeqNum(unsigned int n);

        # Int8 GetEffectiveSearchSpace() const;
        # void SetEffectiveSearchSpace(Int8 eff);
        # void SetEffectiveSearchSpace(const vector<Int8>& eff);

        # int GetDbGeneticCode() const;
        # void SetDbGeneticCode(int gc);

        # const char* GetPHIPattern() const;
        # void SetPHIPattern(const char* pattern, bool is_dna);

        # double GetInclusionThreshold() const
        # void SetInclusionThreshold(double u)

        # int GetPseudoCount() const
        # void SetPseudoCount(int u)

        # bool GetIgnoreMsaMaster() const
        # void SetIgnoreMsaMaster(bool val)
        
        # double GetDomainInclusionThreshold(void) const
        # void SetDomainInclusionThreshold(double th)

        # bool GetUseIndex() const
        # bool GetForceIndex() const
        # bool GetIsOldStyleMBIndex() const
        # bool GetMBIndexLoaded() const
        # const string GetIndexName() const
        # void SetMBIndexLoaded( bool index_loaded = true )
        # void SetUseIndex( 
        #         bool use_index = true, const string & index_name = "", 
        #         bool force_index = false, bool old_style_index = false )

        # void DebugDump(CDebugDumpContext ddc, unsigned int depth) const
        # void DoneDefaults() const
        
        # typedef ncbi::objects::CBlast4_parameters TBlast4Opts
        # TBlast4Opts * GetBlast4AlgoOpts()

        # bool operator==(const CBlastOptions& rhs) const
        # bool operator!=(const CBlastOptions& rhs) const

        # void SetRemoteProgramAndService_Blast3(const string & p, const string & s)
        # void GetRemoteProgramAndService_Blast3(string & p, string & s) const

        # const CBlastOptionsMemento* CreateSnapshot() const
        # void SetDefaultsMode(bool dmode)
        # bool GetDefaultsMode() const