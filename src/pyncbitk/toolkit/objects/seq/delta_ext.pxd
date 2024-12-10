from libcpp cimport bool
from libcpp.list cimport list as cpplist

from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from .delta_seq cimport CDelta_seq

cdef extern from "objects/seq/Delta_ext_.hpp" namespace "ncbi::objects::CDelta_ext_Base":

    ctypedef cpplist[CRef[CDelta_seq]] Tdata


cdef extern from "objects/seq/Delta_ext_.hpp" namespace "ncbi::objects":
    
    cppclass CDelta_ext_Base(CSerialObject):
        CDelta_ext_Base()

        bool IsSet() const
        bool CanGet() const
        void Reset() except +
        
        const Tdata& Get() except +
        Tdata& GetMut "Set"() except +

        # operator const Tdata& void() except +
        # operator Tdata& void() except +


cdef extern from "objects/seq/Delta_ext.hpp" namespace "ncbi::objects":

    cppclass CDelta_ext(CDelta_ext_Base):
        CDelta_ext()

        # CDelta_seq& AddLiteral(TSeqPos len) except +
        # CDelta_seq& AddLiteral(const CTempString& iupac_seq, CSeq_inst::EMol mol) except +
        # CDelta_seq& AddLiteral(const CTempString& iupac_seq, CSeq_inst::EMol mol, bool do_pack) except +

        # void AddAndSplit(const CTempString& src, CSeq_data::E_Choice format, TSeqPos length) except +
        # void AddAndSplit(const CTempString& src, CSeq_data::E_Choice format, TSeqPos length, bool gaps_ok) except + 
        # void AddAndSplit(const CTempString& src, CSeq_data::E_Choice format, TSeqPos length, bool gaps_ok, bool allow_packing) except +

        # CDelta_seq& AddSeqRange(const CSeq_id& id, TSeqPos from_, TSeqPos to) except +
        # CDelta_seq& AddSeqRange(const CSeq_id& id, TSeqPos from_, TSeqPos to, ENa_strand strand) except +
