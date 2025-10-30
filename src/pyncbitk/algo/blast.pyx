# coding: utf-8
# cython: language_level=3

from libcpp cimport bool
from libcpp.cast cimport static_cast
from libcpp.string cimport string
from iostream cimport istream, ostream, filebuf

from ..toolkit.algo.blast.core.blast_options cimport BLAST_DEFAULT_MATRIX, BLAST_GENETIC_CODE
from ..toolkit.algo.blast.api.bl2seq cimport CBl2Seq
from ..toolkit.algo.blast.api.blast_types cimport EProgram, ProgramNameToEnum, TSeqAlignVector, EProgramToTaskName, TSearchMessages
from ..toolkit.algo.blast.api.sseqloc cimport SSeqLoc, TSeqLocVector, CBlastSearchQuery, CBlastQueryVector
from ..toolkit.algo.blast.api.local_blast cimport CLocalBlast
from ..toolkit.algo.blast.api.blast_options cimport CBlastOptions, EAPILocality
from ..toolkit.algo.blast.api.blast_options_handle cimport CBlastOptionsHandle, CBlastOptionsFactory
from ..toolkit.algo.blast.api.blast_nucl_options cimport CBlastNucleotideOptionsHandle
from ..toolkit.algo.blast.api.blast_prot_options cimport CBlastProteinOptionsHandle
from ..toolkit.algo.blast.api.blastx_options cimport CBlastxOptionsHandle
from ..toolkit.algo.blast.api.blast_advprot_options cimport CBlastAdvancedProteinOptionsHandle
from ..toolkit.algo.blast.api.tblastn_options cimport CTBlastnOptionsHandle
from ..toolkit.algo.blast.api.query_data cimport IQueryFactory
from ..toolkit.algo.blast.api.objmgr_query_data cimport CObjMgr_QueryFactory
from ..toolkit.algo.blast.api.objmgrfree_query_data cimport CObjMgrFree_QueryFactory
from ..toolkit.algo.blast.api.local_db_adapter cimport CLocalDbAdapter
from ..toolkit.algo.blast.api.blast_results cimport CSearchResultSet, CSearchResults, size_type as CSearchResults_size_type
from ..toolkit.algo.blast.format.blast_format cimport CBlastFormat
from ..toolkit.algo.blast.blastinput.blast_args cimport EOutputFormat
from ..toolkit.corelib.ncbiobj cimport CConstRef, CRef
from ..toolkit.corelib.ncbistre cimport CNcbiOstream
from ..toolkit.corelib.ncbistr cimport kEmptyStr
from ..toolkit.objects.seq.bioseq cimport CBioseq
from ..toolkit.objects.seqset.bioseq_set cimport CBioseq_set
from ..toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from ..toolkit.objects.seqloc.seq_id cimport CSeq_id
from ..toolkit.objects.seqloc.seq_loc cimport CSeq_loc
from ..toolkit.objmgr.object_manager cimport CObjectManager
from ..toolkit.objmgr.scope cimport CScope
from ..toolkit.objtools.align_format.format_flags cimport kDfltArgNumAlignments, kDfltArgNumDescriptions
from ..toolkit.objtools.readers.fasta cimport CFastaReader
from ..toolkit.serial.serialbase cimport CSerialObject, MSerial_Format_AsnText
from ..toolkit.serial.serialdef cimport ESerialRecursionMode
from ..toolkit.algo.blast.api.uniform_search cimport CSearchDatabase, EMoleculeType
from ..toolkit.objtools.blast.seqdb_reader.seqdb cimport ESeqType

from ..objects.general cimport ObjectId
from ..objects.seqloc cimport SeqLoc
from ..objects.seqid cimport SeqId
from ..objects.seqalign cimport SeqAlign, SeqAlignSet
from ..objects.seq cimport BioSeq
from ..objects.seqset cimport BioSeqSet
from ..objmgr cimport Scope
from ..objtools cimport DatabaseReader

from pystreambuf cimport pywritebuf

from .._utils import peekable, is_iterable


# --- BLAST input --------------------------------------------------------------

cdef class SearchQuery:

    def __init__(self, SeqLoc seqloc not None, Scope scope not None):
        cdef CBlastSearchQuery* _query = new CBlastSearchQuery(seqloc._loc.GetObject(), scope._scope.GetObject())
        self._query.Reset(_query)

    @property
    def seqloc(self):
        cdef CConstRef[CSeq_loc] cref = self._query.GetObject().GetQuerySeqLoc()
        cdef CRef[CSeq_loc] ref = CRef[CSeq_loc](&cref.GetObject())
        return SeqLoc._wrap(ref)

    @property
    def length(self):
        return self._query.GetObject().GetLength()

    @property
    def scope(self):
        cdef Scope scope = Scope.__new__(Scope)
        scope._scope = self._query.GetObject().GetScope()
        return scope


cdef class SearchQueryVector:

    def __init__(self, queries = ()):
        cdef SearchQuery         query
        cdef CBlastQueryVector*  qv    = new CBlastQueryVector()
        for query in queries:
            qv.AddQuery(query._query)
        self._qv.Reset(qv)

    def __len__(self):
        return self._qv.GetObject().Size()


# --- BLAST results ------------------------------------------------------------

cdef class SearchResultsSet:  # TODO: change name to SearchResultsList?

    @staticmethod
    cdef SearchResultsSet _wrap(CRef[CSearchResultSet] ref):
        cdef SearchResultsSet obj = SearchResultsSet.__new__(SearchResultsSet)
        obj._ref = ref
        return obj

    def __len__(self):
        return self._ref.GetNonNullPointer().size()

    def __getitem__(self, ssize_t index):
        cdef ssize_t _index  = index
        cdef ssize_t _length = self._ref.GetNonNullPointer().size()

        if _index < 0:
            _index += _length
        if _index < 0 or _index >= _length:
            raise IndexError(index)

        cdef CSearchResultSet* obj  = self._ref.GetNonNullPointer()
        cdef CSearchResults*   item = &obj[0][<CSearchResults_size_type> _index]
        return SearchResults._wrap(CRef[CSearchResults](item))


cdef class SearchResults:

    @staticmethod
    cdef SearchResults _wrap(CRef[CSearchResults] ref):
        cdef SearchResults obj = SearchResults.__new__(SearchResults)
        obj._ref = ref
        return obj

    @property
    def query_id(self):
        cdef CSeq_id* seq_id = <CSeq_id*> self._ref.GetNonNullPointer().GetSeqId().GetNonNullPointer()
        return SeqId._wrap(CRef[CSeq_id](seq_id))

    @property
    def alignments(self):
        return SeqAlignSet._wrap(self._ref.GetNonNullPointer().GetSeqAlignMut())


# --- BLAST --------------------------------------------------------------------

cdef class Blast:
    """A base command object for running a BLAST search.
    """

    @staticmethod
    def tasks():
        return [ x.decode() for x in CBlastOptionsFactory.GetTasks() ]

    def __init__(
        self,
        *,
        object evalue = None,
        object gapped = None,
        object window_size = None,
        object max_target_sequences = None,
        object xdrop_gap = None,
        object culling_limit = None,
        object percent_identity = None,
    ):
        """__init__(self, *, gapped=None, window_size=None, evalue=None, max_target_sequences=None, culling_limit=None, percent_identity=None)\n--\n
        """
        if self._opt.Empty():
            raise TypeError("Cannot instantiate abstract class Blast")

        if evalue is not None:
            self.evalue = evalue
        if max_target_sequences is not None:
            self.max_target_sequences = max_target_sequences
        if window_size is not None:
            self.window_size = window_size
        if gapped is not None:
            self.gapped = gapped
        if xdrop_gap is not None:
            self.xdrop_gap = xdrop_gap
        if culling_limit is not None:
            self.culling_limit = culling_limit
        if percent_identity is not None:
            self.percent_identity = percent_identity

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.program!r}, window_size={self.window_size})"

    # --- Properties -----------------------------------------------------------

    @property
    def program(self):
        """`str`: The name of the BLAST program.
        """
        cdef EProgram program = self._opt.GetNonNullPointer().GetOptions().GetProgram()
        return EProgramToTaskName(program).decode('ascii')

    @property
    def window_size(self):
        """`int`: The window size for multiple hits.
        """
        return self._opt.GetNonNullPointer().GetWindowSize()

    @window_size.setter
    def window_size(self, int window_size):
        if window_size < 0:
            raise ValueError(f"window_size must be a positive integer, got {window_size!r}")
        self._opt.GetNonNullPointer().SetWindowSize(window_size)

    @property
    def off_diagonal_range(self):
        """`int`: The number of off-diagonals to search for the second hit.
        """
        return self._opt.GetNonNullPointer().GetOffDiagonalRange()

    @property
    def xdrop_gap(self):
        """`float`: X-dropoff value (in bits) for preliminary gapped extensions.
        """
        return self._opt.GetNonNullPointer().GetGapXDropoff()

    @xdrop_gap.setter
    def xdrop_gap(self, double xdrop_gap):
        self._opt.GetNonNullPointer().SetGapXDropoff(xdrop_gap)

    @property
    def xdrop_gap_final(self):
        """`float`: X-dropoff value (in bits) for final gapped alignment.
        """
        return self._opt.GetNonNullPointer().GetGapXDropoffFinal()

    @property
    def evalue(self):
        """`float`: Expectation value (E) threshold for saving hits.
        """
        return self._opt.GetNonNullPointer().GetEvalueThreshold()

    @evalue.setter
    def evalue(self, double evalue):
        if evalue <= 0:
            raise ValueError(f"`evalue` must be greater than zero, got {evalue!r}")
        self._opt.GetNonNullPointer().SetEvalueThreshold(evalue)

    @property
    def percent_identity(self):
        """`float`: Percentage identity threshold for saving hits.
        """
        return self._opt.GetNonNullPointer().GetPercentIdentity()

    @percent_identity.setter
    def percent_identity(self, float percent_identity):
        if percent_identity < 0.0 or percent_identity > 100.0:
            raise ValueError("invalid `percent_identity`: {percent_identity!r}")
        self._opt.GetNonNullPointer().SetPercentIdentity(percent_identity)

    @property
    def coverage_hsp(self):
        """`float`: Query coverage percentage per HSP.
        """
        return self._opt.GetNonNullPointer().GetQueryCovHspPerc()

    @property
    def gapped(self):
        """`bool`: `False` if alignments are performed in ungapped mode only.
        """
        return self._opt.GetNonNullPointer().GetGappedMode()

    @gapped.setter
    def gapped(self, bool gapped):
        self._opt.GetNonNullPointer().SetGappedMode(gapped)

    @property
    def culling_limit(self):
        """`int`: The culling limit for hits.

        If the query range of a hit is enveloped by that of at least this many
        higher-scoring hits, delete the hit.

        """
        return self._opt.GetNonNullPointer().GetCullingLimit()

    @culling_limit.setter
    def culling_limit(self, int culling_limit):
        if culling_limit < 0:
            raise ValueError(f"invalid `culling_limit`: {culling_limit!r}")
        self._opt.GetNonNullPointer().SetCullingLimit(culling_limit)

    @property
    def database_size(self):
        """`int`: The effective length of the database.
        """
        return self._opt.GetNonNullPointer().GetDbLength()

    @property
    def search_space(self):
        """`int`: The effective length of the search space.
        """
        return self._opt.GetNonNullPointer().GetEffectiveSearchSpace()

    @property
    def max_target_sequences(self):
        """`int`: The maximum number of aligned sequences to retain.
        """
        return self._opt.GetNonNullPointer().GetHitlistSize()

    @max_target_sequences.setter
    def max_target_sequences(self, int max_target_sequences):
        if max_target_sequences <= 0:
            raise ValueError(f"`max_target_sequences` must be greater than zero, got {max_target_sequences!r}")
        self._opt.GetNonNullPointer().SetHitlistSize(max_target_sequences)

    # --- Private Methods ------------------------------------------------------

    cdef CBlastOptionsHandle* _options(self) except NULL:
        # TODO: check if this works with anything (as only the CBlastOptions is
        #       accessed) or whether this requires a dedicated implementation
        #       for each CBlastOptionsHandle subclass.
        cdef CBlastNucleotideOptionsHandle* opt = new CBlastNucleotideOptionsHandle(
            self._opt.GetObject().GetOptions().Clone()
        )
        return <CBlastOptionsHandle*> opt

    # --- Public Methods -------------------------------------------------------

    cpdef SearchResultsSet run(
        self,
        BlastQueries queries,
        BlastSubjects subjects,
        bool pairwise = False
    ):
        """Run a BLAST query with the given sequences.

        Arguments:
            queries (`BioSeq`, `BioSeqSet`, `SearchQuery`, `SearchQueryVector`):
                The queries to use on the subject sequences.
            subjects (`BioSeq`, `BioSeqSet`, `SearchQuery`, `SearchQueryVector`,
                `DatabaseReader`): The subjects sequences to search. A BLAST
                database can be given by passing a
                `~pyncbitk.objtools.DatabaseReader` object directly.
            pairwise (`bool`): Set to `True` to run the database search in
                pairwise mode, forcing BLAST to produce one `SearchResults`
                per query, even when no hits were found. *Ignored when*
                ``subjects`` *is a* `DatabaseReader`.

        Returns:
            `~pyncbitk.algo.SearchResultsSet`: The list of search results,
            with one `~pyncbitk.algo.SearchResults` item per query.

        """
        cdef CBlastQueryVector     _queries_loc
        cdef CBlastQueryVector     _subjects_loc
        cdef SearchQuery           query
        cdef CRef[IQueryFactory]   query_factory
        cdef CRef[IQueryFactory]   subject_factory
        cdef CRef[CLocalDbAdapter] db
        cdef CRef[CLocalBlast]     blast

        cdef CRef[CBlastOptionsHandle] opt
        opt.Reset(self._options())

        # prepare queries: a list of `SearchQuery` objects
        if BlastQueries is BioSeq:
            if queries._ref.GetNonNullPointer().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = queries.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            query_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq](queries._ref)))
        elif BlastQueries is BioSeqSet:
            query_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq_set](queries._ref)))
        elif BlastQueries is SearchQuery:
            return self.run(SearchQueryVector((queries,)), subjects)
        elif BlastQueries is SearchQueryVector:
            if queries._qv.GetObject().Empty():
                raise ValueError("Empty query vector")
            query_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(queries._qv.GetObject()))
        else:
            if not is_iterable(queries):
                queries = (queries, )
            for query in queries:
                _queries_loc.AddQuery(query._query)
            query_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_queries_loc))

        # prepare subjects: either a list of `SearchQuery` objects, or a `DatabaseReader`
        if BlastSubjects is DatabaseReader:
            _ty = subjects._ref.GetNonNullPointer().GetSequenceType()
            if _ty == ESeqType.eProtein:
                search_database = new CSearchDatabase(string(b"protein_db"), EMoleculeType.eBlastDbIsProtein)
            elif _ty == ESeqType.eNucleotide:
                search_database = new CSearchDatabase(string(b"nucleotide_db"), EMoleculeType.eBlastDbIsNucleotide)
            else:
                raise ValueError(f"invalid sequence type: {_ty!r}")
            search_database.SetSeqDb(subjects._ref)
            db.Reset(new CLocalDbAdapter(search_database[0]))
        elif BlastSubjects is BioSeq:
            if subjects._ref.GetNonNullPointer().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = subjects.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq](subjects._ref)))
            db.Reset(new CLocalDbAdapter(subject_factory, opt, not pairwise))
        elif BlastSubjects is BioSeqSet:
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq_set](subjects._ref)))
            db.Reset(new CLocalDbAdapter(subject_factory, opt, not pairwise))
        elif BlastSubjects is SearchQuery:
            return self.run(queries, SearchQueryVector((subjects,)))
        elif BlastSubjects is SearchQueryVector:
            if subjects._qv.GetObject().Empty():
                raise ValueError("Empty subjects vector")
            subject_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(subjects._qv.GetObject()))
            db.Reset(new CLocalDbAdapter(subject_factory, opt, not pairwise))
        else:
            if not is_iterable(subjects):
                subjects = (subjects, )
            for query in subjects:
                _subjects_loc.AddQuery(query._query)
            subject_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_subjects_loc))
            db.Reset(new CLocalDbAdapter(subject_factory, opt, not pairwise))

        # prepare the BLAST program
        try:
            with nogil:
                blast.Reset(new CLocalBlast(query_factory, opt, db))
        except RuntimeError as err:
            raise RuntimeError("Failed initializing BLAST") from err
        # if (m_InterruptFnx != NULL) {
        #     m_Blast->SetInterruptCallback(m_InterruptFnx, m_InterruptUserData);
        # }
        # // Set the hitlist size to the total number of subject sequences, to
        # // make sure that no hits are discarded (ported from CBl2Seq::SetupSearch
        # m_OptsHandle.SetHitlistSize((int) m_tSubjects.size());

        with nogil:
            results = blast.GetNonNullPointer().Run()

        # check for warnings or errors
        messages = blast.GetNonNullPointer().GetSearchMessages()
        if messages.HasMessages():
            print(messages.ToString().decode()) # FIXME

        return SearchResultsSet._wrap(results)


cdef class NucleotideBlast(Blast):
    """A base command object for running a nucleotide BLAST search.
    """


cdef class ProteinBlast(Blast):
    """A base command object for running a protein BLAST search.
    """


cdef class MappingBlast(Blast):
    """A base command object for running a mapping BLAST search.
    """


cdef class BlastP(ProteinBlast):
    """A command object for running ``blastn`` searches.
    """

    def __init__(
        self,
        *,
        object word_threshold = None,
        object word_size = None,
        **kwargs,
    ):
        cdef CBlastAdvancedProteinOptionsHandle* handle = new CBlastAdvancedProteinOptionsHandle()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(**kwargs)
        if word_threshold is not None:
            self.word_threshold = word_threshold
        if word_size is not None:
            self.word_size = word_size

    @property
    def word_threshold(self):
        """`float`: The minimum score to record a word in the BLAST lookup table.
        """
        cdef CBlastOptionsHandle* opt = self._opt.GetNonNullPointer()
        cdef CBlastProteinOptionsHandle* popt = <CBlastProteinOptionsHandle*> opt
        return popt.GetWordThreshold()

    @word_threshold.setter
    def word_threshold(self, double word_threshold):
        cdef CBlastOptionsHandle* opt = self._opt.GetNonNullPointer()
        cdef CBlastProteinOptionsHandle* popt = <CBlastProteinOptionsHandle*> opt
        popt.SetWordThreshold(word_threshold)

    @property
    def word_size(self):
        """`int`: The word size for the wordfinder algorithm.
        """
        cdef CBlastOptionsHandle* opt = self._opt.GetNonNullPointer()
        cdef CBlastProteinOptionsHandle* popt = <CBlastProteinOptionsHandle*> opt
        return popt.GetWordThreshold()

    @word_size.setter
    def word_size(self, int word_size):
        cdef CBlastOptionsHandle* opt = self._opt.GetNonNullPointer()
        cdef CBlastProteinOptionsHandle* popt = <CBlastProteinOptionsHandle*> opt
        popt.SetWordSize(word_size)


cdef class BlastN(NucleotideBlast):
    """A command object for running ``blastn`` searches.
    """

    def __init__(
        self,
        *,
        object dust_filtering = None,
        object penalty = None,
        object reward = None,
        **kwargs,
    ):
        cdef CBlastNucleotideOptionsHandle* handle = new CBlastNucleotideOptionsHandle()
        handle.SetTraditionalBlastnDefaults()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(**kwargs)
        if dust_filtering is not None:
            self.dust_filtering = dust_filtering
        if penalty is not None:
            self.penalty = penalty
        if reward is not None:
            self.reward = reward

    @property
    def dust_filtering(self):
        """`bool`: Whether DUST filtering is enabled or not.
        """
        cdef CBlastNucleotideOptionsHandle* handle = <CBlastNucleotideOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetDustFiltering()

    @dust_filtering.setter
    def dust_filtering(self, bint dust_filtering):
        cdef CBlastNucleotideOptionsHandle* handle = <CBlastNucleotideOptionsHandle*> self._opt.GetNonNullPointer()
        handle.SetDustFiltering(dust_filtering)

    @property
    def penalty(self):
        """`int`: The (negative) score penalty for a nucleotide mismatch.
        """
        cdef CBlastNucleotideOptionsHandle* handle = <CBlastNucleotideOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetMismatchPenalty()

    @penalty.setter
    def penalty(self, int penalty):
        cdef CBlastNucleotideOptionsHandle* handle = <CBlastNucleotideOptionsHandle*> self._opt.GetNonNullPointer()
        handle.SetMismatchPenalty(penalty)

    @property
    def reward(self):
        """`int`: The (positive) score reward for a nucleotide match.
        """
        cdef CBlastNucleotideOptionsHandle* handle = <CBlastNucleotideOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetMatchReward()

    @reward.setter
    def reward(self, int reward):
        cdef CBlastNucleotideOptionsHandle* handle = <CBlastNucleotideOptionsHandle*> self._opt.GetNonNullPointer()
        handle.SetMatchReward(reward)


cdef class BlastX(NucleotideBlast):
    """A command object for running ``blastx`` searches.
    """

    def __init__(
        self,
        *,
        int query_genetic_code = 1,
        int max_intron_length = 0,
        **kwargs,
    ):
        cdef CBlastxOptionsHandle* handle = new CBlastxOptionsHandle()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(**kwargs)
        self.query_genetic_code = query_genetic_code
        self.max_intron_length = max_intron_length

    @property
    def max_intron_length(self):
        """`int`: Largest allowed intron in a translated nucleotide sequence.
        """
        cdef CTBlastnOptionsHandle* handle = <CTBlastnOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetLongestIntronLength()

    @max_intron_length.setter
    def max_intron_length(self, int max_intron_length):
        cdef CTBlastnOptionsHandle* handle = <CTBlastnOptionsHandle*> self._opt.GetNonNullPointer()
        if max_intron_length < 0:
            raise ValueError(f"`max_target_sequences` must be a positive integer, got {max_intron_length!r}")
        handle.SetLongestIntronLength(max_intron_length)

    @property
    def query_genetic_code(self):
        """`int`: Genetic code to use for translating the query sequences.
        """
        cdef CBlastxOptionsHandle* handle = <CBlastxOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetQueryGeneticCode()

    @query_genetic_code.setter
    def query_genetic_code(self, int query_genetic_code):
        cdef CBlastxOptionsHandle* handle = <CBlastxOptionsHandle*> self._opt.GetNonNullPointer()
        handle.SetQueryGeneticCode(query_genetic_code)


cdef class TBlastN(ProteinBlast):
    """A command object for running ``tblastn`` searches.
    """

    def __init__(
        self,
        *,
        int database_genetic_code = 1,
        int max_intron_length = 0,
        **kwargs,
    ):
        cdef CTBlastnOptionsHandle* handle = new CTBlastnOptionsHandle()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(**kwargs)
        self.database_genetic_code = database_genetic_code
        self.max_intron_length = max_intron_length

    @property
    def max_intron_length(self):
        """`int`: Largest allowed intron in a translated nucleotide sequence.
        """
        cdef CTBlastnOptionsHandle* handle = <CTBlastnOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetLongestIntronLength()

    @max_intron_length.setter
    def max_intron_length(self, int max_intron_length):
        cdef CTBlastnOptionsHandle* handle = <CTBlastnOptionsHandle*> self._opt.GetNonNullPointer()
        if max_intron_length < 0:
            raise ValueError(f"`max_target_sequences` must be a positive integer, got {max_intron_length!r}")
        handle.SetLongestIntronLength(max_intron_length)

    @property
    def database_genetic_code(self):
        """`int`: Genetic code to use for translating the database sequences.
        """
        cdef CTBlastnOptionsHandle* handle = <CTBlastnOptionsHandle*> self._opt.GetNonNullPointer()
        return handle.GetDbGeneticCode()

    @database_genetic_code.setter
    def database_genetic_code(self, int database_genetic_code):
        cdef CTBlastnOptionsHandle* handle = <CTBlastnOptionsHandle*> self._opt.GetNonNullPointer()
        handle.SetDbGeneticCode(database_genetic_code)


# --- Formatter ----------------------------------------------------------------

cdef class _BlastFormatter:  # WIP

    def __init__(
        self,
        Blast blast,
        object subjects,
        Scope scope,
    ):
        cdef SeqLoc                seqloc
        cdef CRef[CLocalDbAdapter] db
        cdef CBlastFormat*         fmt
        cdef CRef[IQueryFactory]   subject_factory
        cdef CBlastQueryVector     _subjects_loc
        cdef bool                  scan_mode       = False

        self._file = open("/tmp/out.tsv", "wb")
        self._outfile = new ostream(new pywritebuf(self._file))

        # prepare subjects: either a list of `SearchQuery` objects, or a `DatabaseReader`
        if isinstance(subjects, DatabaseReader):
            _ty = (<DatabaseReader> subjects)._ref.GetNonNullPointer().GetSequenceType()
            if _ty == ESeqType.eProtein:
                search_database = new CSearchDatabase(string(), EMoleculeType.eBlastDbIsProtein)
            elif _ty == ESeqType.eNucleotide:
                search_database = new CSearchDatabase(string(), EMoleculeType.eBlastDbIsNucleotide)
            else:
                raise ValueError(f"invalid sequence type: {_ty!r}")
            search_database.SetSeqDb((<DatabaseReader> subjects)._ref)
            db.Reset(new CLocalDbAdapter(search_database[0]))
        elif isinstance(subjects, BioSeq):
            if (<BioSeq> subjects)._ref.GetNonNullPointer().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = subjects.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq]((<BioSeq> subjects)._ref)))
            db.Reset(new CLocalDbAdapter(subject_factory, blast._opt, scan_mode))
        elif isinstance(subjects, BioSeqSet):
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq_set]((<BioSeqSet> subjects)._ref)))
            db.Reset(new CLocalDbAdapter(subject_factory, blast._opt, scan_mode))
        else:
            # if not is_iterable(subjects):
            #     subjects = (subjects, )
            # for seqloc in subjects:
            #     _subjects_loc.push_back(seqloc._seqloc)
            # subject_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_subjects_loc))
            # db.Reset(new CLocalDbAdapter(subject_factory, blast._opt, scan_mode))
            raise NotImplemented


        # make formatter
        self._fmt.Reset(
            new CBlastFormat(
                blast._opt.GetObject().GetOptions(),
                db.GetObject(),
                EOutputFormat.eTabularWithComments,
                True, # believe_query
                self._outfile[0],
                kDfltArgNumDescriptions, # num_summary
                kDfltArgNumAlignments, # num_alignment,
                scope._scope.GetObject(),

                "",
                False, # show_gi,
                False, # is_html
                BLAST_GENETIC_CODE, # qgencode
                BLAST_GENETIC_CODE, # dbgencode
                False, # use_sum_statistics
                False, # is_remote_search
                -1, # dbfilt_algorithm,
                kEmptyStr, # b"qacc sacc pident qcovs qstart qend sstart send", # custom_output_format,
                False, # is_megablast
                False, # is_indexed
                NULL, # ig_opts,
                NULL, # domain_db_adapter
                kEmptyStr, # cmdline
                kEmptyStr, # subjectTag


            )
        )

    def print(self, SearchResults results, SearchQueryVector queries):
        self._fmt.GetObject().PrintOneResultSet(
            results._ref.GetObject(),
            CConstRef[CBlastQueryVector](queries._qv),
        )