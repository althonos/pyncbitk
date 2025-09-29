from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ....corelib.ncbiobj cimport CObject
from ....corelib.ncbitype cimport Uint8
from ....objects.seq.bioseq cimport CBioseq
from ....objects.blastdb.blast_def_line_set cimport CBlast_def_line_set
from ....objmgr.seq_vector cimport CSeqVector
from ..seqdb_reader.seqdbcommon cimport EBlastDbVersion, EOidMaskType

cdef extern from "objtools/blast/seqdb_writer/writedb.hpp" namespace "ncbi::CWriteDB" nogil:

    enum ESeqType:
        eProtein
        eNucleotide

    enum EIndexType:
        eNoIndex
        eSparseIndex
        eFullIndex
        eAddTrace
        eFullWithTrace
        eDefault
        eAddHash

    ctypedef int TIndexType


cdef extern from "objtools/blast/seqdb_writer/writedb.hpp" namespace "ncbi" nogil:

    cppclass CWriteDB(CObject):

        CWriteDB(
            const string & dbname,
            ESeqType seqtype,
            const string & title
        ) except +

        CWriteDB(
            const string& dbname,
            ESeqType seqtype,
            const string& title,
            int itype, # = EIndexType.eDefault,
            bool parse_ids, # = True,
            bool long_ids, # = False,
            bool use_gi_mask, # = False,
            EBlastDbVersion dbver, # = EBlastDbVersion.eBDB_Version4,
            bool limit_defline, # = False,
            Uint8 oid_masks, # = EOidMaskType.fNone,
            bool scan_bioseq_4_cfastareader_usrobj, # = False
        ) except +

        void AddSequence(const CBioseq & bs) except +
        void AddSequence(const CBioseq & bs, CSeqVector & sv) except +
        # void AddSequence(const CBioseq_Handle & bsh)
        # void AddSequence(const CTempString & sequence, const CTempString & ambiguities = "")
        void SetPig(int pig)
        void SetDeflines(const CBlast_def_line_set & deflines) except +

        # int RegisterMaskAlgorithm(EBlast_filter_program program,
        #                       const string & options = string(),
        #                       const string & name = string())
        # int RegisterMaskAlgorithm(const string & id,
        #                       const string & description = string(),
        #                       const string & options = string())
        # void SetMaskData(const CMaskedRangesVector & ranges,
        #              const vector<TGi>         & gis)


        void ListVolumes(vector[string]& vols) except +
        void ListFiles(vector[string]& files) except +

        void Close() except +

        void SetMaxFileSize(Uint8 sz) except +
        void SetMaxVolumeLetters(Uint8 letters) except +

        # @staticmethod
        # CRef[CBlast_def_line_set] ExtractBioseqDeflines(const CBioseq & bs, bool parse_ids=true,
        #                     bool long_ids=false,
        #                     bool scan_bioseq_4_cfastareader_usrobj=false)

        void SetMaskedLetters(const string & masked) except +

        int FindColumn(const string & title) except +
        int CreateUserColumn(const string & title) except +

        # void AddColumnMetaData(int            col_id,
        #                    const string & key,
        #                    const string & value)

        # CBlastDbBlob & SetBlobData(int column_id)