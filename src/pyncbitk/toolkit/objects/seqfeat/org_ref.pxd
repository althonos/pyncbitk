from libcpp cimport bool
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp.vector cimport vector

from ...serial.serialbase cimport CSerialObject
from ...corelib.ncbiobj cimport CRef
from ..general.dbtag cimport CDbtag


cdef extern from "objects/seqfeat/Org_ref_.hpp" namespace "ncbi::objects::COrg_ref_base" nogil:

    ctypedef string TTaxname
    ctypedef string TCommon
    ctypedef list[string] TMod
    ctypedef vector[CRef[CDbtag]] TDb
    ctypedef list[string] TSyn
    # ctypedef COrgName TOrgname


cdef extern from "objects/seqfeat/Org_ref_.hpp" namespace "ncbi::objects" nogil:

    cppclass COrg_ref_Base(CSerialObject):
        COrg_ref_Base()

        # # preferred formal name
        # # optional
        # # typedef string TTaxname
        # #  Check whether the Taxname data member has been assigned a value.
        # bool IsSetTaxname() const
        # # Check whether it is safe or not to call GetTaxname method.
        # bool CanGetTaxname() const
        # void ResetTaxname()
        # const TTaxname& GetTaxname() const
        # void SetTaxname(const TTaxname& value)
        # void SetTaxname(TTaxname&& value)
        # TTaxname& SetTaxname()

        # # common name
        # # optional
        # # typedef string TCommon
        # #  Check whether the Common data member has been assigned a value.
        # bool IsSetCommon() const
        # # Check whether it is safe or not to call GetCommon method.
        # bool CanGetCommon() const
        # void ResetCommon()
        # const TCommon& GetCommon() const
        # void SetCommon(const TCommon& value)
        # void SetCommon(TCommon&& value)
        # TCommon& SetCommon()

        # # unstructured modifiers
        # # optional
        # # typedef list< string > TMod
        # #  Check whether the Mod data member has been assigned a value.
        # bool IsSetMod() const
        # # Check whether it is safe or not to call GetMod method.
        # bool CanGetMod() const
        # void ResetMod()
        # const TMod& GetMod() const
        # TMod& SetMod()

        # # ids in taxonomic or culture dbases
        # # optional
        # # typedef vector< CRef< CDbtag > > TDb
        # #  Check whether the Db data member has been assigned a value.
        # bool IsSetDb() const
        # # Check whether it is safe or not to call GetDb method.
        # bool CanGetDb() const
        # void ResetDb()
        # const TDb& GetDb() const
        # TDb& SetDb()

        # # synonyms for taxname or common
        # # optional
        # # typedef list< string > TSyn
        # #  Check whether the Syn data member has been assigned a value.
        # bool IsSetSyn() const
        # # Check whether it is safe or not to call GetSyn method.
        # bool CanGetSyn() const
        # void ResetSyn()
        # const TSyn& GetSyn() const
        # TSyn& SetSyn()

        # # optional
        # # typedef COrgName TOrgname
        # #  Check whether the Orgname data member has been assigned a value.
        # bool IsSetOrgname() const
        # # Check whether it is safe or not to call GetOrgname method.
        # bool CanGetOrgname() const
        # void ResetOrgname()
        # const TOrgname& GetOrgname() const
        # void SetOrgname(TOrgname& value)
        # TOrgname& SetOrgname()

        # # Reset the whole object
        # virtual void Reset()



cdef extern from "objects/seqfeat/Org_ref.hpp" namespace "ncbi::objects::COrg_ref" nogil:

    enum EOrgref_part:
        eOrgref_nothing
        eOrgref_taxname
        eOrgref_common 
        eOrgref_mod    
        eOrgref_db
        eOrgref_db_taxid
        eOrgref_db_all  
        eOrgref_syn
        eOrgref_orgname
        eOrgref_on_name
        eOrgref_on_attr
        eOrgref_on_attr_spec
        eOrgref_on_attr_nofwd
        eOrgref_on_attr_uncult
        eOrgref_on_attr_all
        eOrgref_on_mod    
        eOrgref_on_mod_nom
        eOrgref_on_mod_oldname
        eOrgref_on_mod_tm
        eOrgref_on_mod_all
        eOrgref_on_lin
        eOrgref_on_gc
        eOrgref_on_mgc
        eOrgref_on_pgc
        eOrgref_on_div
        eOrgref_on_all
        eOrgref_all
        eOrgref_all_but_syn
        eOrgref_all_but_spec
        eOrgref_default

    ctypedef unsigned int fOrgref_parts


cdef extern from "objects/seqfeat/Org_ref.hpp" namespace "ncbi::objects" nogil:

    cppclass COrg_ref(COrg_ref_Base):
        COrg_ref()

        # # Appends a label to "label" based on content
        # void GetLabel(string* label) const

        # # Returns NCBI Taxonomy database id (AKA tax id)
        # # if the latter is found in Org_ref otherwise returns 0
        # TTaxId GetTaxId() const
        # # Sets tax id into Org_ref contents.
        # # Returns old value of tax id or 0 if it was not found
        # TTaxId SetTaxId( TTaxId tax_id )

        # # shortcut access to selected OrgName methods
        # bool IsSetLineage() const
        # const string& GetLineage() const
        
        # bool IsSetGcode() const
        # int GetGcode() const
        
        # bool IsSetMgcode() const
        # int GetMgcode() const
        
        # bool IsSetPgcode() const
        # int GetPgcode() const
        
        # bool IsSetDivision() const
        # const string& GetDivision() const
        
        # bool IsSetOrgMod() const

        # bool IsVarietyValid(const string& variety) const
        # bool HasValidVariety() const
        # bool IsSubspeciesValid(const string& subspecies) const

        # CRef<COrg_ref> MakeCommon(const COrg_ref& other) const
        
        # static CConstRef<COrg_ref> TableLookup(const string& taxname)
        # static const vector<string>& GetTaxnameList()
        # bool UpdateFromTable()

        # void CleanForGenBank()

        # void FilterOutParts( fOrgref_parts to_remain )
