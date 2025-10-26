from libcpp cimport bool
from libcpp.string cimport string

from .object_id cimport CObject_id
from ...corelib.ncbimisc cimport TTaxId
from ...serial.serialbase cimport CSerialObject

cdef extern from "objects/general/Dbtag_.hpp" namespace "ncbi::objects::CDbtag_Base" nogil:

    ctypedef string TDb
    ctypedef CObject_id TTag


cdef extern from "objects/general/Dbtag_.hpp" namespace "ncbi::objects" nogil:

    cppclass CDbtag_Base(CSerialObject):
        CDbtag_Base()
        void Reset()

        bool IsSetDb() const
        bool CanGetDb() const
        void ResetDb()
        const TDb& GetDb() except +
        void SetDb(const TDb& value) except +
        # void SetDb(TDb&& value) except +
        TDb& GetDbMut "SetDb"() except +

        bool IsSetTag() const
        bool CanGetTag() const
        void ResetTag()
        const TTag& GetTag() except +
        void SetTag(TTag& value) except +
        TTag& GetTagMut "SetTag"() except +


cdef extern from "objects/general/Dbtag.hpp" namespace "ncbi::objects::CDbtag" nogil:

    enum EDbtagType:
        eDbtagType_bad
        eDbtagType_AFTOL
        eDbtagType_APHIDBASE
        eDbtagType_ASAP
        eDbtagType_ATCC
        eDbtagType_ATCC_dna
        eDbtagType_ATCC_in_host
        eDbtagType_AceView_WormGenes
        eDbtagType_AntWeb
        eDbtagType_ApiDB
        eDbtagType_ApiDB_CryptoDB
        eDbtagType_ApiDB_PlasmoDB
        eDbtagType_ApiDB_ToxoDB
        eDbtagType_BB
        eDbtagType_BDGP_EST
        eDbtagType_BDGP_INS
        eDbtagType_BEETLEBASE
        eDbtagType_BGD
        eDbtagType_BoLD
        eDbtagType_CCDS
        eDbtagType_CDD
        eDbtagType_CGNC
        eDbtagType_CK
        eDbtagType_COG
        eDbtagType_CloneID
        eDbtagType_CollecTF
        eDbtagType_DDBJ
        eDbtagType_ECOCYC
        eDbtagType_EMBL
        eDbtagType_ENSEMBL
        eDbtagType_ESTLIB
        eDbtagType_EcoGene
        eDbtagType_FANTOM_DB
        eDbtagType_FBOL
        eDbtagType_FLYBASE
        eDbtagType_Fungorum
        eDbtagType_GABI
        eDbtagType_GDB
        eDbtagType_GEO
        eDbtagType_GI
        eDbtagType_GO
        eDbtagType_GOA
        eDbtagType_GRIN
        eDbtagType_GeneDB
        eDbtagType_GeneID
        eDbtagType_GrainGenes
        eDbtagType_Greengenes
        eDbtagType_HGNC
        eDbtagType_HMP
        eDbtagType_HOMD
        eDbtagType_HPM
        eDbtagType_HPRD
        eDbtagType_HSSP
        eDbtagType_H_InvDB
        eDbtagType_IFO
        eDbtagType_IMGT_GENEDB
        eDbtagType_IMGT_HLA
        eDbtagType_IMGT_LIGM
        eDbtagType_IRD
        eDbtagType_ISD
        eDbtagType_ISFinder
        eDbtagType_InterimID
        eDbtagType_Interpro
        eDbtagType_IntrepidBio
        eDbtagType_JCM
        eDbtagType_JGIDB
        eDbtagType_LRG
        eDbtagType_LocusID
        eDbtagType_MGI
        eDbtagType_MIM
        eDbtagType_MaizeGDB
        eDbtagType_MycoBank
        eDbtagType_NMPDR
        eDbtagType_NRESTdb
        eDbtagType_NextDB
        eDbtagType_OrthoMCL
        eDbtagType_Osa1
        eDbtagType_PBR
        eDbtagType_PBmice
        eDbtagType_PDB
        eDbtagType_PFAM
        eDbtagType_PGN
        eDbtagType_PIR
        eDbtagType_PSEUDO
        eDbtagType_Pathema
        eDbtagType_Phytozome
        eDbtagType_PomBase
        eDbtagType_PseudoCap
        eDbtagType_RAP_DB
        eDbtagType_RATMAP
        eDbtagType_RBGE_garden
        eDbtagType_RBGE_herbarium
        eDbtagType_REBASE
        eDbtagType_RFAM
        eDbtagType_RGD
        eDbtagType_RZPD
        eDbtagType_RiceGenes
        eDbtagType_SEED
        eDbtagType_SGD
        eDbtagType_SGN
        eDbtagType_SK_FST
        eDbtagType_SRPDB
        eDbtagType_SoyBase
        eDbtagType_SubtiList
        eDbtagType_TAIR
        eDbtagType_TIGRFAM
        eDbtagType_UNILIB
        eDbtagType_UNITE
        eDbtagType_UniGene
        eDbtagType_UniProt_SwissProt
        eDbtagType_UniProt_TrEMBL
        eDbtagType_UniSTS
        eDbtagType_VBASE2
        eDbtagType_VBRC
        eDbtagType_VectorBase
        eDbtagType_Vega
        eDbtagType_WorfDB
        eDbtagType_WormBase
        eDbtagType_Xenbase
        eDbtagType_ZFIN
        eDbtagType_axeldb
        eDbtagType_dbClone
        eDbtagType_dbCloneLib
        eDbtagType_dbEST
        eDbtagType_dbProbe
        eDbtagType_dbSNP
        eDbtagType_dbSTS
        eDbtagType_dictyBase
        eDbtagType_miRBase
        eDbtagType_niaEST
        eDbtagType_taxon
        eDbtagType_MGD
        eDbtagType_PID
        eDbtagType_BEEBASE
        eDbtagType_NASONIABASE
        eDbtagType_BioProject
        eDbtagType_IKMC
        eDbtagType_ViPR
        eDbtagType_PubChem
        eDbtagType_SRA
        eDbtagType_Trace
        eDbtagType_RefSeq
        eDbtagType_EnsemblGenomes
        eDbtagType_EnsemblGenomes_Gn
        eDbtagType_EnsemblGenomes_Tr
        eDbtagType_TubercuList
        eDbtagType_MedGen
        eDbtagType_CGD
        eDbtagType_Assembly
        eDbtagType_GenBank
        eDbtagType_BioSample
        eDbtagType_ISHAM_ITS
        eDbtagType_ERIC
        eDbtagType_I5KNAL
        eDbtagType_VISTA
        eDbtagType_BEI
        eDbtagType_Araport
        eDbtagType_VGNC
        eDbtagType_RNAcentral
        eDbtagType_PeptideAtlas
        eDbtagType_EPDnew
        eDbtagType_Ensembl
        eDbtagType_PseudoCAP
        eDbtagType_MarpolBase
        eDbtagType_dbVar
        eDbtagType_EnsemblRapid
        eDbtagType_AllianceGenome
        eDbtagType_EchinoBase
        eDbtagType_AmoebaDB
        eDbtagType_CryptoDB
        eDbtagType_FungiDB
        eDbtagType_GiardiaDB
        eDbtagType_MicrosporidiaDB
        eDbtagType_PiroplasmaDB
        eDbtagType_PlasmoDB
        eDbtagType_ToxoDB
        eDbtagType_TrichDB
        eDbtagType_TriTrypDB
        eDbtagType_VEuPathDB
        eDbtagType_NCBIOrtholog

    enum EDbtagGroup:
        fNone
        fGenBank
        fRefSeq
        fSrc   
        fProbe

    ctypedef int TDbtagGroup

    enum EIsRefseq:
        eIsRefseq_No
        eIsRefseq_Yes
    
    enum EIsSource:
        eIsSource_No
        eIsSource_Yes

    enum EIsEstOrGss:
        eIsEstOrGss_No
        eIsEstOrGss_Yes


cdef extern from "objects/general/Dbtag.hpp" namespace "ncbi::objects" nogil:

    cppclass CDbtag(CDbtag_Base):
        CDbtag()

        int Compare(const CDbtag& dbt2) const
        bool Match(const CDbtag& dbt2) const
        bool SetAsMatchingTo(const CDbtag& dbt2)
        void GetLabel(string* label) const

        # bool IsApproved(EIsRefseq refseq = eIsRefseq_No, EIsSource is_source = eIsSource_No, EIsEstOrGss is_est_or_gss = eIsEstOrGss_No ) const;
        # bool IsApproved(TDbtagGroup group) const

        bool IsSkippable() const
        EDbtagType GetType() const
        bool GetDBFlags (bool& is_refseq, bool& is_src, string& correct_caps) const
        TDbtagGroup GetDBFlags (string& correct_caps) const
        void InvalidateType()

        string GetUrl() const
        string GetUrl(TTaxId taxid) const
        string GetUrl(const string &taxname) const
        # string GetUrl(
        #     const string & genus,
        #     const string & species,
        #     const string & subspecies = kEmptyStr
        # ) const

