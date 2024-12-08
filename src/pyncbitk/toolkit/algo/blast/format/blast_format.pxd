from libcpp cimport bool
from libcpp.string cimport string

from ....corelib.ncbistre cimport CNcbiOstream
from ....corelib.ncbistr cimport kEmptyStr
from ....corelib.ncbiobj cimport CObject, CRef, CConstRef
from ....objmgr.scope cimport CScope
from ..api.blast_options cimport CBlastOptions
from ..api.blast_results cimport CSearchResults
from ..api.local_db_adapter cimport CLocalDbAdapter
from ..api.sseqloc cimport CBlastQueryVector
from ..blastinput.blast_args cimport EOutputFormat
from ..igblast.igblast cimport CIgBlastOptions


cdef extern from "algo/blast/format/blast_format.hpp" namespace "ncbi::CBlastFormat" nogil:

    cppclass SClone:
        SClone()

        string na
        string chain_type
        string aa
        string v_gene
        string d_gene
        string j_gene
        string c_gene
        string seqid
        double identity
        string productive
    
    enum DisplayOption:
        eDescriptions,
        eAlignments,
        eMetadata,
        eDescriptionsWithTemplates

    const int kFormatLineLength
    const int kMinTaxFormatLineLength


cdef extern from "algo/blast/format/blast_format.hpp" namespace "ncbi" nogil:

    cppclass CBlastFormat(CObject):
        CBlastFormat(
            const CBlastOptions& opts, 
            CLocalDbAdapter& db_adapter,
            EOutputFormat format_type, 
            bool believe_query, 
            CNcbiOstream& outfile,
            int num_summary, 
            int num_alignments,
            CScope& scope,

            const char *matrix_name, # = BLAST_DEFAULT_MATRIX,
            bool show_gi, # = false, 
            bool is_html, # = false, 
            int qgencode, # = BLAST_GENETIC_CODE,
            int dbgencode, # = BLAST_GENETIC_CODE,
            bool use_sum_statistics, # = false,
            bool is_remote_search, # = false,
            int dbfilt_algorithm, # = -1,
            const string& custom_output_format, # = kEmptyStr,
            bool is_megablast, # = false,
            bool is_indexed, # = false,
            const CIgBlastOptions* ig_opts,# = NULL,
            const CLocalDbAdapter* domain_db_adapter,# = NULL,
            const string & cmdline,# = kEmptyStr,
            const string& subjectTag,# = kEmptyStr

        ) except +

        # CBlastFormat(
        #     const CBlastOptions& opts, 
        #     blast::CLocalDbAdapter& db_adapter,
        #     blast::CFormattingArgs::EOutputFormat format_type, 
        #     bool believe_query, 
        #     CNcbiOstream& outfile,
        #     int num_summary, 
        #     int num_alignments,
        #     CScope & scope,
        #     const char *matrix_name = BLAST_DEFAULT_MATRIX,
        #     bool show_gi = false, 
        #     bool is_html = false, 
        #     int qgencode = BLAST_GENETIC_CODE,
        #     int dbgencode = BLAST_GENETIC_CODE,
        #     bool use_sum_statistics = false,
        #     bool is_remote_search = false,
        #     int dbfilt_algorithm = -1,
        #     const string& custom_output_format = kEmptyStr,
        #     bool is_megablast = false,
        #     bool is_indexed = false,
        #     const blast::CIgBlastOptions * ig_opts = NULL,
        #     const blast::CLocalDbAdapter* domain_db_adapter = NULL,
        #     const string & cmdline = kEmptyStr,
		#     const string& subjectTag = kEmptyStr
        # ) except +

        # CBlastFormat(
        #     const CBlastOptions& opts, 
        #     const vector< CBlastFormatUtil::SDbInfo >& dbinfo_list,
        #     blast::CFormattingArgs::EOutputFormat format_type, 
        #     bool believe_query, CNcbiOstream& outfile,
        #     int num_summary, 
        #     int num_alignments,
        #     CScope & scope,
        #     bool show_gi = false, 
        #     bool is_html = false, 
        #     bool is_remote_search = false,
        #     const string& custom_output_format = kEmptyStr,
        #     bool is_vdb = false,
        #     const string & cmdline = kEmptyStr
        # ) except +

        # void PrintProlog() except ++
        # Int8 GetDbTotalLength() except ++

        void PrintOneResultSet(const CSearchResults& results, CConstRef[CBlastQueryVector]& queries)