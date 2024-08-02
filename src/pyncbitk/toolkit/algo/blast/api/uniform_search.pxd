from libcpp cimport bool
from libcpp.string cimport string

from ....corelib.ncbiobj cimport CObject, CRef
from ....objtools.blast.seqdb_reader.seqdb cimport CSeqDB


cdef extern from "algo/blast/api/uniform_search.hpp" namespace "ncbi::blast::CSearchDatabase" nogil:

    # ctypedef vector[TGi] TGiList

    enum EMoleculeType:
        eBlastDbIsProtein
        eBlastDbIsNucleotide


cdef extern from "algo/blast/api/uniform_search.hpp" namespace "ncbi::blast" nogil:

    cppclass CSearchDatabase(CObject):
        CSearchDatabase(const string& dbname, EMoleculeType mol_type)
        CSearchDatabase(const string& dbname, EMoleculeType mol_type, const string& entrez_query)

        string GetDatabaseName() const
        void SetDatabaseName(const string& dbname)

        void SetMoleculeType(EMoleculeType mol_type)
        EMoleculeType GetMoleculeType() const

        bool IsProtein() const

        void SetEntrezQueryLimitation(const string& entrez_query)
        string GetEntrezQueryLimitation() const

        # void SetGiList(CSeqDBGiList * gilist)
        # const CRef<CSeqDBGiList>& GetGiList() const
        # const TGiList GetGiListLimitation() const
        # void SetNegativeGiList(CSeqDBGiList * gilist)
        # const CRef<CSeqDBGiList>& GetNegativeGiList() const
        # const TGiList GetNegativeGiListLimitation() const

        # void SetFilteringAlgorithm(int filt_algorithm_id, ESubjectMaskingType mask_type);
        # void SetFilteringAlgorithm(const string &filt_algorithm, ESubjectMaskingType mask_type);
        # int GetFilteringAlgorithm() const;
        # string GetFilteringAlgorithmKey() const;
        # ESubjectMaskingType GetMaskType() const;

        void SetSeqDb(CRef[CSeqDB] seqdb)
        CRef[CSeqDB] GetSeqDb() const