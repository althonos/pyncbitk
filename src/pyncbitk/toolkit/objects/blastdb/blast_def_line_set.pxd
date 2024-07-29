from libcpp cimport bool
from libcpp.list cimport list

from ...corelib.ncbiobj cimport CRef
from ...serial.serialbase cimport CSerialObject
from ...corelib.ncbimisc cimport TGi
from .blast_def_line cimport CBlast_def_line


cdef extern from "objects/blastdb/Blast_def_line_set_.hpp" namespace "ncbi::objects::CBlast_def_line_set_Base" nogil:

    ctypedef list[CRef[CBlast_def_line]] TData


cdef extern from "objects/blastdb/Blast_def_line_set_.hpp" namespace "ncbi::objects" nogil:

    cppclass CBlast_def_line_set_Base(CSerialObject):
        CBlast_def_line_set_Base()


cdef extern from "objects/blastdb/Blast_def_line_set.hpp" namespace "ncbi::objects" nogil:

    cppclass CBlast_def_line_set(CBlast_def_line_set_Base):
        CBlast_def_line_set()

        void SortBySeqIdRank(bool is_protein)
        void SortBySeqIdRank(bool is_protein, bool useBlastRank)
        void PutTargetGiFirst(TGi gi)
        void RemoveGIs()
