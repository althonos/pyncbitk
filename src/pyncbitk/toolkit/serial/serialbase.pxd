from libcpp cimport bool
from libcpp.string cimport string

from ..corelib.ncbiobj cimport CObject, CRef
from ..corelib.ncbistre cimport CNcbiIstream, CNcbiOstream
from .serialdef cimport ESerialRecursionMode, ESerialDataFormat, TSerial_AsnText_Flags, TSerial_Json_Flags, TSerial_Xml_Flags

cdef extern from "serial/serialbase.hpp" namespace "ncbi" nogil:

    enum EResetVariant:
        eDoResetVariant
        eDoNotResetVariant

    cppclass CSerialObject(CObject):
        CSerialObject()

        void Assign(const CSerialObject& source) except +
        void Assign(const CSerialObject& source, ESerialRecursionMode how) except +

        bool Equals(const CSerialObject& object) except +
        bool Equals(const CSerialObject& object, ESerialRecursionMode how) except +

        CNcbiOstream& operator<< (CNcbiOstream& str, const CSerialObject&    obj)
        CNcbiIstream& operator>> (CNcbiIstream& str, CSerialObject&          obj)


    cppclass CAliasBase[TPrim]:
        CAliasBase()
        CAliasBase(const TPrim& value)
        
        const TPrim& Get() const
        TPrim& GetMut "Set" ()
        void Set(const TPrim& value)

    cppclass CStringAliasBase[TString](CAliasBase[TString]):
        CStringAliasBase()
        CStringAliasBase(const TString& value)


    cppclass MSerial_Flags:
        pass

    CNcbiOstream& operator<< (CNcbiOstream& io, const MSerial_Flags& obj) except +
    CNcbiIstream& operator>> (CNcbiIstream& io, const MSerial_Flags& obj) except +
    

    ctypedef unsigned int TSerial_Format_Flags
    cppclass MSerial_Format(MSerial_Flags):
        MSerial_Format(ESerialDataFormat fmt)
        MSerial_Format(ESerialDataFormat fmt, TSerial_Format_Flags flags)

    cppclass MSerial_Format_AsnText(MSerial_Format):
        MSerial_Format_AsnText()
        MSerial_Format& operator()(TSerial_AsnText_Flags flags)

    cppclass MSerial_Format_AsnBinary(MSerial_Format):
        MSerial_Format_AsnBinary()

    cppclass MSerial_Format_Xml(MSerial_Format):
        MSerial_Format_Xml()
        MSerial_Format& operator()(TSerial_Xml_Flags flags)

    cppclass MSerial_Format_Json(MSerial_Format):
        MSerial_Format_Json()
        MSerial_Format& operator()(TSerial_Json_Flags flags)


    const char* operator>>(const char* s, CSerialObject& obj) except +
    string operator>>(const string& s, CSerialObject& obj) except +
    string& operator<<(string& s, const CSerialObject& obj) except +