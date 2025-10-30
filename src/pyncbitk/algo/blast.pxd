# coding: utf-8
# cython: language_level=3

from libcpp cimport bool

from ..toolkit.algo.blast.api.sseqloc cimport CBlastSearchQuery, CBlastQueryVector
from ..toolkit.algo.blast.api.blast_options_handle cimport CBlastOptionsHandle
from ..toolkit.algo.blast.api.blast_results cimport CSearchResultSet, CSearchResults
from ..toolkit.algo.blast.format.blast_format cimport CBlastFormat
from ..toolkit.corelib.ncbiobj cimport CRef
from ..toolkit.corelib.ncbistre cimport CNcbiOstream

from ..objects.seq cimport BioSeq
from ..objects.seqset cimport BioSeqSet
from ..objtools cimport DatabaseReader

# --- BLAST input --------------------------------------------------------------

cdef class SearchQuery:
    cdef CRef[CBlastSearchQuery] _query


cdef class SearchQueryVector:
    cdef CRef[CBlastQueryVector] _qv


ctypedef fused BlastQueries:
    BioSeq
    BioSeqSet
    SearchQuery
    SearchQueryVector
    object

ctypedef fused BlastSubjects:
    BioSeq
    BioSeqSet
    DatabaseReader
    SearchQuery
    SearchQueryVector
    object

# --- BLAST results ------------------------------------------------------------

cdef class SearchResultsSet:  # TODO: change name to SearchResultsList?
    cdef CRef[CSearchResultSet] _ref

    @staticmethod
    cdef SearchResultsSet _wrap(CRef[CSearchResultSet] ref)


cdef class SearchResults:
    cdef CRef[CSearchResults] _ref

    @staticmethod
    cdef SearchResults _wrap(CRef[CSearchResults] ref)


# --- BLAST --------------------------------------------------------------------

cdef class Blast:
    """A base command object for running a BLAST search.
    """
    cdef CRef[CBlastOptionsHandle] _opt

    cdef CBlastOptionsHandle* _options(self) except NULL

    cpdef SearchResultsSet run(
        self,
        BlastQueries queries,
        BlastSubjects subjects,
        bool pairwise = *,
    )


cdef class NucleotideBlast(Blast):
    pass

cdef class ProteinBlast(Blast):
    pass

cdef class MappingBlast(Blast):
    pass

cdef class BlastP(ProteinBlast):
    pass
    

cdef class BlastN(NucleotideBlast):
    pass


cdef class BlastX(NucleotideBlast):
    pass


cdef class TBlastN(ProteinBlast):
    pass


# --- Formatter ----------------------------------------------------------------

cdef class _BlastFormatter:  # WIP
    cdef CRef[CBlastFormat] _fmt
    cdef object             _file
    cdef CNcbiOstream*      _outfile
