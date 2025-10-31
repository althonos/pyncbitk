from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.ncbitype cimport Uint8
from ...corelib.ncbidiag cimport EDiagSev
from ...corelib.ncbistre cimport CNcbiOstream


cdef extern from "objtools/logging/message.hpp" namespace "ncbi::objects" nogil:

    cppclass IObjtoolsMessage:
        IObjtoolsMessage *Clone() except +
        void Write(CNcbiOstream& out) except +
        void Dump(CNcbiOstream& out) except +
        void WriteAsXML(CNcbiOstream& out) except +
        void DumpAsXML(CNcbiOstream& out) except +
        string GetText() except +
        EDiagSev GetSeverity() except +
        int GetCode() except +
        int GetSubCode() except +

    cppclass CObjtoolsMessage(IObjtoolsMessage):
        CObjtoolsMessage(const string& text, EDiagSev severity)
        CObjtoolsMessage *Clone() except +
        # string Compose() const
        void Write(CNcbiOstream& out) except +
        void Dump(CNcbiOstream& out) except +
        void WriteAsXML(CNcbiOstream& out) except +
        void DumpAsXML(CNcbiOstream& out) except +
        string GetText() except +
        EDiagSev GetSeverity() except +
        int GetCode() except +
        int GetSubCode() except+
