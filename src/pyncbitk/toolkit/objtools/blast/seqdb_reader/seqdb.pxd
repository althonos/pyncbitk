from libcpp cimport bool
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from ....corelib.ncbiobj cimport CObject, CRef
from ....corelib.ncbimisc cimport TGi
from ....corelib.ncbitype cimport Int8, Uint8
from ....objects.seq.bioseq cimport CBioseq
from ....objects.seqloc.seq_id cimport CSeq_id
from .seqdbcommon cimport EBlastDbVersion

cdef extern from "objtools/blast/seqdb_reader/seqdb.hpp" namespace "ncbi::CSeqDB" nogil:

    enum EOidListType:
        eOidList
        eOidRange

    enum ESeqType:
        eProtein
        eNucleotide
        eUnknown

    enum ESummaryType:
        eUnfilteredAll
        eFilteredAll
        eFilteredRange

    enum EMmapFileTypes:
        eMmap_IndexFile
        eMmap_SequenceFile

    enum EMmapStrategies:
        eMmap_Normal
        eMmap_Sequential
        eMmap_WillNeed

    ctypedef int TOID
    ctypedef int TPIG
    ctypedef TGi TGI

    # struct TOffsetPair:
    #     TSeqPos first
    #     TSeqPos second

    const char* kBlastDbDateFormat


cdef extern from "objtools/blast/seqdb_reader/seqdb.hpp" namespace "ncbi" nogil:

    cppclass CSeqDBIter:
        CSeqDBIter operator++()

        int GetOID()
        const char* GetData()
        int GetLength()

        CSeqDBIter(const CSeqDBIter &)
        CSeqDBIter& operator=(const CSeqDBIter &)

        bool operator!() const 

    cppclass CSeqDB(CObject):
        CSeqDB(const string& dbname, ESeqType seqtype) except +
        # CSeqDB(const string& dbname, ESeqType seqtype, CSeqDBGiList* gilist = NULL, bool use_atlas_lock = True) except +
        # CSeqDB(const string& dbname, ESeqType seqtype, CSeqDBNegativeList* nlist) except +
        # CSeqDB(const string& dbname, ESeqType seqtype, CSeqDBGiList* gilist, CSeqDBNegativeList* nlist) except +
        # CSeqDB(const string& dbname, ESeqType seqtype, int oid_begin, int oid_end, CSeqDBGiList* gilist, CSeqDBNegativeList* nlist) except +
        # CSeqDB(const string& dbname, ESeqType seqtype, CSeqDBIdSet ids) except +
        # CSeqDB(const vector[string]& dbs, ESeqType seqtype, CSeqDBGiList* gilist = NULL) except +
        # CSeqDB(const string& dbname, ESeqType seqtype, int oid_begin, int oid_end, bool use_mmap, CSeqDBGiList* gi_list = NULl) except +
        # CSeqDB(const vector[string]& dbname, ESeqType seqtype, int oid_begin, int oid_end, bool use_mmap, CSeqDBGiList* gi_list = NULl) except +

        @staticmethod
        string GenerateSearchPath() except +

        # int GetSeqLength(int oid) const
        # TGi GetSeqGI(int oid) const
        # int GetSeqLengthApprox(int oid)

        # CRef[CBlast_def_line_set] GetHdr(int oid) const
        # void GetLeafTaxIDs(int oid, map[TGi, set[TTaxId]]& gi_to_taxid_set, bool persist = False)
        # void GetLeafTaxIDs(int oid, vector[TTaxId]& taxids, bool persist = False)
        # void GetTaxIDs(int oid, map[TGi, TTaxId]& gi_to_taxid, bool persist = False)
        # void GetTaxIDs(int oid, vector[TTaxId]& taxids, bool persist = False)
        # void GetAllTaxIds(int oid, set[TTaxId]& taxids)
        CRef[CBioseq] GetBioseq(int oid) except +
        CRef[CBioseq] GetBioseq(int oid, TGi target_gi) except +
        CRef[CBioseq] GetBioseq(int oid, TGi target_gi, const CSeq_id* target_seq_id) except +
        CRef[CBioseq] GetBioseqNoData(int oid, TGi target_gi = ZERO_GI, const CSeq_id* target_set_id = NULL) except +

        # int GetSequence(int oid, const char** buffer) const
        # int GetAmbigSeq(int oid, const char** buffer, int nucl_code) const
        # int GetAmbigSeq(int oid, const char** buffer, int nucl_code, int begin_offset, int end_offset) const
        # int GetAmbigSeqAlloc(int oid, char** buffer, int nucl_code, ESeqDBAllocType strategy, TSequenceRanges* masks = NULL) const
        # int GetAmbigPartialSeq(int oid, char** buffer, int nucl_code, ESeqDBAllocType strategy, TSequenceRanges* partial_ranges, TSequenceRanges* masks = NULL) const
        # void RetSequence(const char** buffer) const
        # void RetSequence(const char** buffer) const
        list[CRef[CSeq_id]] GetSeqIDs(int oid) const
        # void GetGis(int oid, vector[TGi] & gis, bool append = false) const

        ESeqType GetSequenceType() const
        string GetTitle() const
        string GetDate() const

        # @staticmethod
        # CTime GetDate(const string& dbname, ESeqType seqtype)

        int GetNumSeqs() const
        int GetNumSeqStats() const
        int GetNumOIDs() const

        Uint8 GetTotalLength() const
        Uint8 GetExactTotalLength() const
        Uint8 GetTotalLengthStats() const
        Uint8 GetVolumeLength() const

        # void GetTotals(ESummaryType   sumtype,
        #            int          * oid_count,
        #            Uint8        * total_length,
        #            bool           use_approx = true) const

        int GetMaxLength() const
        int GetMinLength() const

        CSeqDBIter Begin() const

        bool CheckOrFindOID(int& next_oid) const

        # EOidListType
        # GetNextOIDChunk(int         & begin_chunk,       // out
        #                 int         & end_chunk,         // out
        #                 int         oid_size,            // in
        #                 vector<int> & oid_list,          // out
        #                 int         * oid_state = NULL); // in+out

        void ResetInternalChunkBookmark()

        const string &GetDBNameList() const
        # const CSeqDBGiList* GetGiList() const
        # CSeqDBIdSet GetIdSet() const

        # bool PigToOid(int pig, int & oid) const
        # bool OidToPig(int oid, int & pig) const
        # bool TiToOid(Int8 ti, int & oid) const
        # bool OidToGi(int oid, TGi & gi) const
        # bool GiToOid(TGi gi, int & oid) const
        # bool GiToOidwFilterCheck(TGi gi, int & oid) const
        # bool GiToPig(TGi gi, int & pig) const
        # bool PigToGi(int pig, TGi & gi) const
        # void AccessionToOids(const string & acc, vector<int> & oids) const
        # void AccessionsToOids(const vector<string>& accs, vector<blastdb::TOid>& oids) const
        void SeqidToOids(const CSeq_id& seqid, vector[int]& oids) const
        bool SeqidToOid(const CSeq_id& seqid, int& oid) const

        # int GetOidAtOffset(int first_seq, Uint8 residue) const
        # CRef[CBioseq] GiToBioseq(TGi gi) const
        # CRef[CBioseq] PigToBioseq(int pig) const
        CRef[CBioseq] SeqidToBioseq(const CSeq_id & seqid) const

        # @staticmethod
        # FindVolumePaths(const string   & dbname,
        #                 ESeqType         seqtype,
        #                 vector<string> & paths,
        #                 vector<string> * alias_paths = NULL,
        #                 bool             recursive = true,
        #                 bool             expand_links = true);

        # void FindVolumePaths(vector<string> & paths, bool recursive=true) const

        # void SetIterationRange(int oid_begin, int oid_end)
        # void GetAliasFileValues(TAliasFileValues & afv)

        # @staticmethod
        # void GetTaxInfo(TTaxId taxid, SSeqDBTaxInfo & info)

        # CRef<CSeq_data> GetSeqData(int oid, TSeqPos begin, TSeqPos end) except +
        # void GetSequenceAsString(int oid, CSeqUtil::ECoding coding, string& output, TSeqRange range = TSeqRange())
        # void GetSequenceAsString(int oid, string& output, TSeqRange range = TSeqRange()) const

        void ListColumns(vector[string] & titles)
        int GetColumnId(const string & title)
        const map[string,string] & GetColumnMetaData(int column_id)
        const string & GetColumnValue(int column_id, const string & key)
        const map[string,string] & GetColumnMetaData(int column_id, const string & volname)
        # void GetColumnBlob(int col_id, int oid, CBlastDbBlob & blob)
        void GetAvailableMaskAlgorithms(vector[int] & algorithms)
        int GetMaskAlgorithmId(const string &algo_name) const
        string GetAvailableMaskAlgorithmDescriptions()
        vector[int] ValidateMaskAlgorithms(const vector[int]& algorithm_ids)

        # void GetMaskAlgorithmDetails(int algorithm_id, 
        #     objects::EBlast_filter_program & program,
        #     string& program_name,
        #     string& algo_opts
        # )
        # void GetMaskAlgorithmDetails(
        #     int algorithm_id,
        #     string& program,
        #     string& program_name,
        #     string& algo_opts
        # )
        # void GetMaskData(
        #     int oid,
        #     const vector<int> & algo_ids,
        #     TSequenceRanges   & ranges
        # )
        # void GetMaskData(
        #     int              oid,
        #     int              algo_id,
        #     TSequenceRanges &ranges
        # )

        # void SetOffsetRanges(int oid, const TRangeList & offset_ranges, bool append_ranges, bool cache_data)
        void RemoveOffsetRanges(int oid)
        void FlushOffsetRangeCache()

        void SetNumberOfThreads(int num_threads, bool force_mt = false)
        Int8 GetDiskUsage() const
        void SetVolsMemBit(int mbit)
        # void DebugDump(CDebugDumpContext ddc, unsigned int depth) const
        EBlastDbVersion GetBlastDbVersion() const

        # void TaxIdsToOids(set<TTaxId>& tax_ids, vector<blastdb::TOid>& rv) const
        # void GetDBTaxIds(set<TTaxId> & tax_ids) const
        # void GetTaxIdsForOids(const vector<blastdb::TOid> & oids, set<TTaxId> & tax_ids) const

        # CRef<CBlast_db_metadata> GetDBMetaData(string user_path = kEmptyStr)

        # void GetTaxIdsForAccession(const string & accs, vector<TTaxId> & taxids)
        # void GetTaxIdsForSeqId(const CSeq_id & seq_id, vector<TTaxId> & taxids)