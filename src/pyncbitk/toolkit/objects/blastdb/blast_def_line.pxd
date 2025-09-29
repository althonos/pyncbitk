from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.set cimport set as cppset
from libcpp.string cimport string

from ...corelib.ncbiobj cimport CRef
from ...corelib.ncbimisc cimport TTaxId as _TTaxId
from ...serial.serialbase cimport CSerialObject
from ..seqloc.seq_id cimport CSeq_id


cdef extern from "objects/blastdb/Blast_def_line_.hpp" namespace "ncbi::objects::CBlast_def_line_Base" nogil:

    ctypedef string TTitle
    ctypedef cpplist[CRef[CSeq_id]] TSeqid
    ctypedef _TTaxId TTaxid
    ctypedef cpplist[int] TMemberships
    ctypedef cpplist[int] TLinks
    ctypedef cpplist[int] TOther_info


cdef extern from "objects/blastdb/Blast_def_line_.hpp" namespace "ncbi::objects" nogil:

    cppclass CBlast_def_line_Base(CSerialObject):
        CBlast_def_line_Base()

        bool IsSetTitle() const
        bool CanGetTitle() const
        void ResetTitle()
        const TTitle& GetTitle() const
        void SetTitle(const TTitle& value)
        void SetTitle(TTitle&& value)
        TTitle& GetTitleMut "SetTitle" ()

        bool IsSetSeqid() const
        bool CanGetSeqid() const
        void ResetSeqid()
        const TSeqid& GetSeqid() const
        TSeqid& GetSeqidMut "SetSeqid" ()

        bool IsSetTaxid() const
        bool CanGetTaxid() const
        void ResetTaxid()
        TTaxid GetTaxid() const
        void SetTaxid(TTaxid value)
        TTaxid& GetTaxidMut "SetTaxid" ()

        bool IsSetMemberships() const
        bool CanGetMemberships() const
        void ResetMemberships()
        const TMemberships& GetMemberships() const
        TMemberships& GetMembershipsMut "SetMemberships" ()

        bool IsSetLinks() const
        bool CanGetLinks() const
        void ResetLinks()
        const TLinks& GetLinks() const
        TLinks& GetLinksMut "SetLinks" ()

        bool IsSetOther_info() const
        bool CanGetOther_info() const
        void ResetOther_info()
        const TOther_info& GetOther_info() const
        TOther_info& GetOther_infoMut "SetOther_info" ()




cdef extern from "objects/blastdb/Blast_def_line.hpp" namespace "ncbi::objects::CBlast_def_line" nogil:

    # ctypedef set[TTaxId] TTaxIds
    pass


cdef extern from "objects/blastdb/Blast_def_line.hpp" namespace "ncbi::objects" nogil:

    cppclass CBlast_def_line(CBlast_def_line_Base):
        CBlast_def_line()

        # TTaxIds GetLeafTaxIds() const
        # void SetLeafTaxIds(const TTaxIds& t)
