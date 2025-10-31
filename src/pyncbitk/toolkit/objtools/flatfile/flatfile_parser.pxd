from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.ncbiobj cimport CRef
from ...corelib.ncbistre cimport CNcbiIstream
from ...serial.serialbase cimport CSerialObject
from ..logging.listener cimport IObjtoolsListener
from .flatfile_parse_info cimport Parser

cdef extern from "objtools/flatfile/flatfile_parser.hpp" namespace "ncbi" nogil:

    cppclass CFlatFileParser:
        CFlatFileParser(IObjtoolsListener* pMessageListener) except +
        #CRef[CSerialObject] Parse(Parser& parseInfo) except +
        CRef[CSerialObject] Parse(Parser& parseInfo, CNcbiIstream& istr) except +
        #CRef[CSerialObject] Parse(Parser& parseInfo, const string& filename) except +
        # bool Parse(Parser& parseInfo, CObjectOStream& objOstr) except +
