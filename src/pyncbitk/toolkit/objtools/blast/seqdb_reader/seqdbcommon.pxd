from libcpp cimport bool
from libcpp.string cimport string

from ....corelib.ncbitype cimport Int4, Uint4, Uint8


cdef extern from "objtools/blast/seqdb_reader/seqdbcommon.hpp" namespace "ncbi::blastdb" nogil:

    ctypedef Int4 TOid


cdef extern from "objtools/blast/seqdb_reader/seqdbcommon.hpp" namespace "ncbi" nogil:

    enum EBlastDbVersion:
        eBDB_Version4
        eBDB_Version5

    const string kSeqDBGroupAliasFileName
    const int kSeqDBNuclNcbiNA8
    const int kSeqDBNuclBlastNA8

    const TOid kSeqDBEntryNotFound
    const TOid kSeqDBEntryDuplicate

    enum ESeqDBAllocType:
        eAtlas
        eMalloc
        eNew

    ctypedef Uint8 TTi
    ctypedef Uint4 TPig

    cppclass SBlastSeqIdListInfo:
        SBlastSeqIdListInfo()   

        bool is_v4
        Uint8 file_size
        Uint8 num_ids
        string title
        string create_date
        Uint8 db_vol_length
        string db_create_date
        string db_vol_names

    enum EOidMaskType:
        fNone
        fExcludeModel

    