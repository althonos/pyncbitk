from libcpp cimport bool
from libcpp.string cimport string

from ...serial.serialbase cimport CSerialObject



cdef extern from "objects/blastdb/Blast_def_line_.hpp" namespace "ncbi::objects::CBlast_def_line_Base" nogil:

    ctypedef string TTitle
    # typedef list< CRef< CSeq_id > > TSeqid
    # typedef NCBI_NS_NCBI::TTaxId TTaxid
    # typedef list< int > TMemberships
    # typedef list< int > TLinks
    # typedef list< int > TOther_info

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



cdef extern from "objects/blastdb/Blast_def_line.hpp" namespace "ncbi::objects" nogil:

    cppclass CBlast_def_line(CBlast_def_line_Base):
        pass
