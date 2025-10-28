from libcpp cimport bool

from ...corelib.tempstr cimport CTempString
from ...serial.serialbase cimport CSerialObject
from ...serial.serializable cimport CSerializable
from ..general.object_id cimport CObject_id
from ..general.dbtag cimport CDbtag
from .textseq_id cimport CTextseq_id

cdef extern from "objects/seqloc/Seq_id_.hpp" namespace "ncbi::objects::CSeq_id_Base" nogil:

    enum E_Choice:
        e_not_set
        e_Local
        e_Gibbsq
        e_Gibbmt
        e_Giim
        e_Genbank
        e_Embl
        e_Pir
        e_Swissprot
        e_Patent
        e_Other
        e_General
        e_Gi
        e_Ddbj
        e_Prf
        e_Pdb
        e_Tpg
        e_Tpe
        e_Tpd
        e_Gpipe
        e_Named_annot_track

    ctypedef CObject_id TLocal
    # typedef int TGibbsq;
    # typedef int TGibbmt;
    # typedef CGiimport_id TGiim;
    ctypedef CTextseq_id TGenbank
    ctypedef CTextseq_id TEmbl
    ctypedef CTextseq_id TPir
    ctypedef CTextseq_id TSwissprot
    # typedef CPatent_seq_id TPatent;
    ctypedef CTextseq_id TOther
    ctypedef CDbtag TGeneral
    # typedef NCBI_NS_NCBI::TGi TGi;
    ctypedef CTextseq_id TDdbj
    ctypedef CTextseq_id TPrf
    # ctypedef CPDB_seq_id TPdb
    ctypedef CTextseq_id TTpg
    ctypedef CTextseq_id TTpe
    ctypedef CTextseq_id TTpd
    ctypedef CTextseq_id TGpipe
    ctypedef CTextseq_id TNamed_annot_track


cdef extern from "objects/seqloc/Seq_id_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_id_Base(CSerialObject):
        CSeq_id_Base()

        void Reset()
        void ResetSelection()

        E_Choice Which() const

        void Select(E_Choice index)
        # void Select(E_Choice index, EResetVariant reset)
        # void Select(E_Choice index, EResetVariant reset, CObjectMemoryPool* pool)

        bool IsLocal() const
        const TLocal& GetLocal() except +
        TLocal& GetLocalMut "SetLocal"() except +
        void SetLocal(TLocal& value) except +

        # bool IsGibbsq(void) const;
        # TGibbsq GetGibbsq(void) const;
        # TGibbsq& SetGibbsq(void);
        # void SetGibbsq(TGibbsq value);

        # bool IsGibbmt(void) const;
        # TGibbmt GetGibbmt(void) const;
        # TGibbmt& SetGibbmt(void);
        # void SetGibbmt(TGibbmt value);

        # bool IsGiim(void) const;
        # const TGiim& GetGiim(void) const;
        # TGiim& SetGiim(void);
        # void SetGiim(TGiim& value);

        bool IsGenbank() const
        const TGenbank& GetGenbank() const
        TGenbank& GetGenbankMut "SetGenbank" ()
        void SetGenbank(TGenbank& value)

        # bool IsEmbl(void) const;
        # const TEmbl& GetEmbl(void) const;
        # TEmbl& SetEmbl(void);
        # void SetEmbl(TEmbl& value);

        # bool IsPir(void) const;
        # const TPir& GetPir(void) const;
        # TPir& SetPir(void);
        # void SetPir(TPir& value);

        # bool IsSwissprot(void) const;
        # const TSwissprot& GetSwissprot(void) const;
        # TSwissprot& SetSwissprot(void);
        # void SetSwissprot(TSwissprot& value);

        # bool IsPatent(void) const;
        # const TPatent& GetPatent(void) const;
        # TPatent& SetPatent(void);
        # void SetPatent(TPatent& value);

        bool IsOther() const
        const TOther& GetOther() except +
        TOther& GetOtherMut "SetOther"() except +
        void SetOther(TOther& value) except +

        bool IsGeneral() const
        const TGeneral& GetGeneral() except +
        TGeneral& GetGeneralMut "SetGeneral" () except +
        void SetGeneral(TGeneral& value) except +

        # bool IsGi(void) const;
        # TGi GetGi(void) const;
        # TGi& SetGi(void);
        # void SetGi(TGi value);

        # bool IsDdbj(void) const;
        # const TDdbj& GetDdbj(void) const;
        # TDdbj& SetDdbj(void);
        # void SetDdbj(TDdbj& value);

        # bool IsPrf(void) const;
        # const TPrf& GetPrf(void) const;
        # TPrf& SetPrf(void);
        # void SetPrf(TPrf& value);

        # bool IsPdb(void) const;
        # const TPdb& GetPdb(void) const;
        # TPdb& SetPdb(void);
        # void SetPdb(TPdb& value);

        # bool IsTpg(void) const;
        # const TTpg& GetTpg(void) const;
        # TTpg& SetTpg(void);
        # void SetTpg(TTpg& value);

        # bool IsTpe(void) const;
        # const TTpe& GetTpe(void) const;
        # TTpe& SetTpe(void);
        # void SetTpe(TTpe& value);

        # bool IsTpd(void) const;
        # const TTpd& GetTpd(void) const;
        # TTpd& SetTpd(void);
        # void SetTpd(TTpd& value);

        # bool IsGpipe(void) const;
        # const TGpipe& GetGpipe(void) const;
        # TGpipe& SetGpipe(void);
        # void SetGpipe(TGpipe& value);

        # bool IsNamed_annot_track(void) const;
        # const TNamed_annot_track& GetNamed_annot_track(void) const;
        # TNamed_annot_track& SetNamed_annot_track(void);
        # void SetNamed_annot_track(TNamed_annot_track& value);


cdef extern from "objects/seqloc/Seq_id_.hpp" namespace "ncbi::objects::CSeq_id" nogil:

    enum EParseFlags:
        fParse_PartialOK
        fParse_RawText  
        fParse_RawGI    
        fParse_AnyRaw   
        fParse_ValidLocal
        fParse_AnyLocal
        fParse_NoFASTA 
        fParse_FallbackOK
        fParse_Default

    ctypedef int TParseFlags

    enum E_SIC:
        e_error
        e_DIFF
        e_NO
        e_YES

cdef extern from "objects/seqloc/Seq_id.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_id(CSeq_id_Base, CSerializable):
        CSeq_id()
        CSeq_id(const CTempString& the_id) except +
        CSeq_id(const CTempString& the_id, TParseFlags flags) except +

        E_SIC Compare(const CSeq_id& sid2) const
