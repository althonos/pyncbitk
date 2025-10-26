from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CSerialObject
from .seqdesc cimport CSeqdesc


cdef extern from "objects/seq/Seqdesc_.hpp" namespace "ncbi::objects::CSeqdesc_Base" nogil:

    enum E_Choice:
        e_not_set
        e_Mol_type
        e_Modif
        e_Method
        e_Name
        e_Title
        e_Org
        e_Comment
        e_Num
        e_Maploc
        e_Pir
        e_Genbank
        e_Pub
        e_Region
        e_User
        e_Sp
        e_Dbxref
        e_Embl
        e_Create_date
        e_Update_date
        e_Prf
        e_Pdb
        e_Het
        e_Source
        e_Molinfo
        e_Modelev

    # ctypedef EGIBB_mol TMol_type
    # ctypedef list< EGIBB_mod > TModif
    # ctypedef EGIBB_method TMethod
    ctypedef string TName
    ctypedef string TTitle
    # ctypedef COrg_ref TOrg
    # ctypedef string TComment
    # ctypedef CNumbering TNum
    # ctypedef CDbtag TMaploc
    # ctypedef CPIR_block TPir
    # ctypedef CGB_block TGenbank
    # ctypedef CPubdesc TPub
    ctypedef string TRegion
    # ctypedef CUser_object TUser
    # ctypedef CSP_block TSp
    # ctypedef CDbtag TDbxref
    # ctypedef CEMBL_block TEmbl
    # ctypedef CDate TCreate_date
    # ctypedef CDate TUpdate_date
    # ctypedef CPRF_block TPrf
    # ctypedef CPDB_block TPdb
    # ctypedef CHeterogen THet
    # ctypedef CBioSource TSource
    # ctypedef CMolInfo TMolinfo
    # ctypedef CModelEvidenceSupport TModelev



cdef extern from "objects/seq/Seqdesc_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeqdesc_Base(CSerialObject):
        CSeqdesc_Base()

        void Reset()
        void ResetSelection()
        E_Choice Which() const
        void CheckSelected(E_Choice index) except +
        void ThrowInvalidSelection(E_Choice) except +
        @staticmethod
        string SelectionName(E_Choice index) noexcept
        void Select(E_Choice index)
        # void Select(E_Choice index, 
        #             EResetVariant reset = eDoResetVariant)
        # void Select(E_Choice index,
        #             EResetVariant reset,
        #             CObjectMemoryPool* pool)
     
cdef extern from "objects/seq/Seqdesc.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeqdesc(CSeqdesc_Base):
        CSeqdesc()

        # bool IsMol_type() const
        # TMol_type GetMol_type() except +
        # TMol_type& SetMol_type() except +
        # void SetMol_type(TMol_type value) except +

        # bool IsModif(void) const;
        # const TModif& GetModif(void) except +
        # TModif& SetModif(void) except +

        # bool IsMethod(void) const;
        # TMethod GetMethod(void) except +
        # TMethod& SetMethod(void) except +
        # void SetMethod(TMethod value) except +

        bool IsName() const
        const TName& GetName() except +
        TName& GetNameMut "SetName" () except +
        void SetName(const TName& value) except +

        bool IsTitle() const
        const TTitle& GetTitle() except +
        TTitle& GetTitleMut "SetTitle" () except +
        void SetTitle(const TTitle& value) except +
    
        # bool IsOrg(void) const;
        # const TOrg& GetOrg(void) const;
        # TOrg& SetOrg(void);
        # void SetOrg(TOrg& value);

        # bool IsComment(void) const;
        # const TComment& GetComment(void) const;
        # TComment& SetComment(void);
        # void SetComment(const TComment& value);

        # bool IsNum(void) const;
        # const TNum& GetNum(void) const;
        # TNum& SetNum(void);
        # void SetNum(TNum& value);

        # bool IsMaploc(void) const;
        # const TMaploc& GetMaploc(void) const;
        # TMaploc& SetMaploc(void);
        # void SetMaploc(TMaploc& value);

        # bool IsPir(void) const;
        # const TPir& GetPir(void) const;
        # TPir& SetPir(void);
        # void SetPir(TPir& value);

        # bool IsGenbank(void) const;
        # const TGenbank& GetGenbank(void) const;
        # TGenbank& SetGenbank(void);
        # void SetGenbank(TGenbank& value);

        # bool IsPub(void) const;
        # const TPub& GetPub(void) const;
        # TPub& SetPub(void);
        # void SetPub(TPub& value);

        bool IsRegion() const
        const TRegion& GetRegion() except +
        TRegion& GetRegionMut "SetRegion" () except +
        void SetRegion(const TRegion& value) except +

        # bool IsUser(void) const;
        # const TUser& GetUser(void) const;
        # TUser& SetUser(void);
        # void SetUser(TUser& value);

        # bool IsSp(void) const;
        # const TSp& GetSp(void) const;
        # TSp& SetSp(void);
        # void SetSp(TSp& value);

        # bool IsDbxref(void) const;
        # const TDbxref& GetDbxref(void) const;
        # TDbxref& SetDbxref(void);
        # void SetDbxref(TDbxref& value);

        # bool IsEmbl(void) const;
        # const TEmbl& GetEmbl(void) const;
        # TEmbl& SetEmbl(void);
        # void SetEmbl(TEmbl& value);

        # bool IsCreate_date(void) const;
        # const TCreate_date& GetCreate_date(void) const;
        # TCreate_date& SetCreate_date(void);
        # void SetCreate_date(TCreate_date& value);

        # bool IsUpdate_date(void) const;
        # const TUpdate_date& GetUpdate_date(void) const;
        # TUpdate_date& SetUpdate_date(void);
        # void SetUpdate_date(TUpdate_date& value);

        # bool IsPrf(void) const;
        # const TPrf& GetPrf(void) const;
        # TPrf& SetPrf(void);
        # void SetPrf(TPrf& value);

        # bool IsPdb(void) const;
        # const TPdb& GetPdb(void) const;
        # TPdb& SetPdb(void);
        # void SetPdb(TPdb& value);

        # bool IsHet(void) const;
        # const THet& GetHet(void) const;
        # THet& SetHet(void);
        # void SetHet(const THet& value);

        # bool IsSource(void) const;
        # const TSource& GetSource(void) const;
        # TSource& SetSource(void);
        # void SetSource(TSource& value);

        # bool IsMolinfo(void) const;
        # const TMolinfo& GetMolinfo(void) const;
        # TMolinfo& SetMolinfo(void);
        # void SetMolinfo(TMolinfo& value);

        # bool IsModelev(void) const;
        # const TModelev& GetModelev(void) const;
        # TModelev& SetModelev(void);
        # void SetModelev(TModelev& value);