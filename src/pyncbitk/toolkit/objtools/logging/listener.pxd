from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.ncbitype cimport Uint8
from ...corelib.ncbidiag cimport EDiagSev
from ...corelib.ncbistre cimport CNcbiOstream
from .message cimport IObjtoolsMessage


cdef extern from "objtools/logging/listener.hpp" namespace "ncbi::objects" nogil:

    cppclass IObjtoolsListener:
        pass

    cppclass CObjtoolsListener(IObjtoolsListener):
        CObjtoolsListener()
        bool PutMessage(const IObjtoolsMessage& message) except +
        void PutProgress(const string& message, const Uint8 iNumDone, const Uint8 iNumTotal) except +
        const IObjtoolsMessage& GetMessage(size_t index) except +
        size_t Count() const
        void ClearAll() except +
        size_t LevelCount(EDiagSev severity) const
        void Dump(CNcbiOstream& ostr) except +
        void DumpAsXML(CNcbiOstream& ostr) except +
        void SetProgressOstream(CNcbiOstream* pProgressOstream) except +

        # class CConstIterator : public TBaseIterator {
        # public:
        #     using value_type = TBaseIterator::value_type::element_type;
        #     using pointer = value_type*;
        #     using reference = value_type&;

        #     CConstIterator(const TBaseIterator& base_it) : TBaseIterator(base_it) {}

        #     reference operator*() const { return *(this->TBaseIterator::operator*()); }
        #     pointer operator->() const { return this->TBaseIterator::operator*().get(); }
        # };    

        # using TConstIterator = CConstIterator;
        # TConstIterator begin(void) const;
        # TConstIterator end(void) const;