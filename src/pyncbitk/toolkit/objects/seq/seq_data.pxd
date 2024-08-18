from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ...corelib.ncbimisc cimport TSeqPos
from ...serial.serialbase cimport CSerialObject
from .iupacna cimport CIUPACna
from .iupacaa cimport CIUPACaa
from .ncbi2na cimport CNCBI2na
from .ncbi4na cimport CNCBI4na
from .ncbi8na cimport CNCBI8na
from .ncbieaa cimport CNCBIeaa

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
    ctypedef CIUPACaa TIupacaa
    ctypedef CNCBI2na TNcbi2na
    ctypedef CNCBI4na TNcbi4na
    ctypedef CNCBI8na TNcbi8na
    # ctypedef CNCBIpna TNcbipna
    # ctypedef CNCBI8aa TNcbi8aa
    ctypedef CNCBIeaa TNcbieaa
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
        #             EResetVariant reset = eDoResetVariant)
        # void Select(E_Choice index,
        #             EResetVariant reset,
        #             CObjectMemoryPool* pool)
        
        bool IsIupacna() const
        const TIupacna& GetIupacna() const
        TIupacna& GetIupacnaMut "SetIupacna" ()
        void SetIupacna(const TIupacna& value)

        bool IsIupacaa() const
        const TIupacaa& GetIupacaa() const
        TIupacaa& GetIupacaaMut "SetIupacaa"()
        void SetIupacaa(const TIupacaa& value)

        bool IsNcbi2na() const
        const TNcbi2na& GetNcbi2na() const
        TNcbi2na& GetNcbi2naMut "SetNcbi2na"()

        bool IsNcbi4na() const
        const TNcbi4na& GetNcbi4na() const
        TNcbi4na& GetNcbi4naMut "SetNcbi4na"()

        bool IsNcbi8na() const
        const TNcbi8na& GetNcbi8na() const
        TNcbi8na& GetNcbi8naMut "SetNcbi8na"()

        # bool IsNcbipna() const
        # const TNcbipna& GetNcbipna() const
        # TNcbipna& GetNcbipnaMut "SetNcbipna"()

        # bool IsNcbi8aa() const
        # const TNcbi8aa& GetNcbi8aa() const
        # TNcbi8aa& GetNcbi8aaMut "SetNcbi8aa"()

        bool IsNcbieaa() const
        const TNcbieaa& GetNcbieaa() const
        TNcbieaa& GetNcbieaaMut "SetNcbieaa"()
        void SetNcbieaa(const TNcbieaa& value)

        # bool IsNcbipaa() const
        # const TNcbipaa& GetNcbipaa() const
        # TNcbipaa& GetNcbipaaMut "SetNcbipaa"()

        # bool IsNcbistdaa() const
        # const TNcbistdaa& GetNcbistdaa() const
        # TNcbistdaa& GetNcbistdaaMut "SetNcbistdaa"()

        # bool IsGap() const
        # const TGap& GetGap() const
        # TGap& SetGap()
        # void SetGap(TGap& value)


cdef extern from "objects/seq/Seq_data.hpp" namespace "ncbi::objects" nogil:

    cppclass CSeq_data(CSeq_data_Base):
        CSeq_data()
        CSeq_data(const string& value, E_Choice index)
        CSeq_data(const vector[char]& value, E_Choice index)
