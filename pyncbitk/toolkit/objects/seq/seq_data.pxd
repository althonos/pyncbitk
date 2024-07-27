from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ...corelib.ncbimisc cimport TSeqPos
from ...serial.serialbase cimport CSerialObject
from .iupacna cimport CIUPACna

cdef extern from "objects/seq/Seq_data_.hpp" namespace "ncbi::objects::CSeq_data_Base" nogil:
    
    enum E_Choice:
        e_not_set
        e_Iupacna
        e_Iupacaa
        e_Ncbi2na
        e_Ncbi4na
        e_Ncbi8na
        e_Ncbipna
        e_Ncbi8aa
        e_Ncbieaa
        e_Ncbipaa
        e_Ncbistdaa
        e_Gap

    ctypedef CIUPACna TIupacna
    # ctypedef CIUPACaa TIupacaa
    # ctypedef CNCBI2na TNcbi2na
    # ctypedef CNCBI4na TNcbi4na
    # ctypedef CNCBI8na TNcbi8na
    # ctypedef CNCBIpna TNcbipna
    # ctypedef CNCBI8aa TNcbi8aa
    # ctypedef CNCBIeaa TNcbieaa
    # ctypedef CNCBIpaa TNcbipaa
    # ctypedef CNCBIstdaa TNcbistdaa
    # ctypedef CSeq_gap TGap

cdef extern from "objects/seq/Seq_data_.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_data_Base(CSerialObject):
        CSeq_data_Base()

        void Reset()
        void ResetSelection()
        E_Choice Which() const
        void CheckSelected(E_Choice index) except +
        void ThrowInvalidSelection(E_Choice index) except +
        @staticmethod
        string SelectionName(E_Choice index) noexcept
        void Select(E_Choice index) except +
        # void Select(E_Choice index, 
        #             EResetVariant reset = eDoResetVariant);
        # void Select(E_Choice index,
        #             EResetVariant reset,
        #             CObjectMemoryPool* pool);
        
        bool IsIupacna() const
        const TIupacna& GetIupacna() const
        TIupacna& GetIupacnaRw "SetIupacna" ()
        void SetIupacna(const TIupacna& value)


cdef extern from "objects/seq/Seq_data.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_data(CSeq_data_Base):
        CSeq_data()
        CSeq_data(const string& value, E_Choice index)
        CSeq_data(const vector[char]& value, E_Choice index)
