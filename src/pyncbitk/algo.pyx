# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.string cimport string

from .toolkit.algo.blast.api.bl2seq cimport CBl2Seq
from .toolkit.algo.blast.api.blast_types cimport EProgram, ProgramNameToEnum, TSeqAlignVector, EProgramToTaskName, TSearchMessages
from .toolkit.algo.blast.api.sseqloc cimport SSeqLoc, TSeqLocVector
from .toolkit.algo.blast.api.local_blast cimport CLocalBlast
from .toolkit.algo.blast.api.blast_options_handle cimport CBlastOptionsHandle, CBlastOptionsFactory
from .toolkit.algo.blast.api.blast_nucl_options cimport CBlastNucleotideOptionsHandle
from .toolkit.algo.blast.api.blast_prot_options cimport CBlastProteinOptionsHandle
from .toolkit.algo.blast.api.blastx_options cimport CBlastxOptionsHandle
from .toolkit.algo.blast.api.blast_advprot_options cimport CBlastAdvancedProteinOptionsHandle
from .toolkit.algo.blast.api.tblastn_options cimport CTBlastnOptionsHandle
from .toolkit.algo.blast.api.query_data cimport IQueryFactory
from .toolkit.algo.blast.api.objmgr_query_data cimport CObjMgr_QueryFactory
from .toolkit.algo.blast.api.objmgrfree_query_data cimport CObjMgrFree_QueryFactory
from .toolkit.algo.blast.api.local_db_adapter cimport CLocalDbAdapter
from .toolkit.algo.blast.api.blast_results cimport CSearchResultSet, CSearchResults, size_type as CSearchResults_size_type
from .toolkit.corelib.ncbiobj cimport CConstRef, CRef
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objects.seqset.bioseq_set cimport CBioseq_set
from .toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from .toolkit.objects.seqloc.seq_id cimport CSeq_id
from .toolkit.objmgr.object_manager cimport CObjectManager
from .toolkit.objmgr.scope cimport CScope
from .toolkit.objtools.readers.fasta cimport CFastaReader
from .toolkit.serial.serialbase cimport CSerialObject, MSerial_Format_AsnText
from .toolkit.serial.serialdef cimport ESerialRecursionMode
from .toolkit.algo.blast.api.uniform_search cimport CSearchDatabase, EMoleculeType
from .toolkit.objtools.blast.seqdb_reader.seqdb cimport ESeqType

from .objects.general cimport ObjectId
from .objects.seqloc cimport SeqLoc, SeqId
from .objects.seqalign cimport SeqAlign, SeqAlignSet
from .objects.seq cimport BioSeq
from .objects.seqset cimport BioSeqSet
from .objmgr cimport Scope
from .objtools cimport DatabaseReader

import os
from ._utils import peekable, is_iterable


# --- BLAST input --------------------------------------------------------------

# cdef class BlastDbLoc:
#     cdef CRef[CSearchDatabase] _ref

#     def __init__(self, str name not None, protein=False):
#         cdef bytes            _name = name.encode() # FIXME: os.fsencode?
#         cdef CSearchDatabase* _db   = new CSearchDatabase(_name, EMoleculeType.eBlastDbIsNucleotide)
#         self._ref.Reset(_db)


cdef class BlastSeqLoc:
    cdef SSeqLoc _seqloc

    def __init__(self, SeqLoc loc, Scope scope, SeqLoc mask = None):
        self._seqloc = SSeqLoc(loc._loc.GetObject(), scope._scope.GetObject())


ctypedef fused BlastQueries:
    BioSeq
    BioSeqSet
    object

ctypedef fused BlastSubjects:
    BioSeq
    BioSeqSet
    DatabaseReader
    object

# --- BLAST results ------------------------------------------------------------

cdef class SearchResultsSet:
    cdef CRef[CSearchResultSet] _ref

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
    cdef CRef[CSearchResults] _ref

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
    cdef CRef[CBlastOptionsHandle] _opt

    @staticmethod
    def tasks():
        return [ x.decode() for x in CBlastOptionsFactory.GetTasks() ]

    def __init__(
        self,
        *,
        bool gapped = True,
        int window_size = 40,
        double evalue = 10.0,
        int max_target_sequences = 500,
    ):
        if self._opt.Empty():
            raise TypeError("Cannot instantiate abstract class Blast")
        self.window_size = window_size
        self.evalue = evalue
        self.max_target_sequences = max_target_sequences
        self.gapped = gapped

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

    # --- Public Methods -------------------------------------------------------

    cpdef SearchResultsSet run(
        self,
        BlastQueries queries,
        BlastSubjects subjects,
        bool scan_mode = False
    ):
        """Run a BLAST query with the given sequences.

        Arguments:
            queries (`BioSeq` or `BioSeqSet`): The queries to use on the
                subject sequences.
            subjects (`BioSeq`, `BioSeqSet`, or `DatabaseReader`): The
                subjects sequences to search. A BLAST database can be given
                by passing a `~pyncbitk.objtools.DatabaseReader` objects
                directly.
            scan_mode (`bool`): Set to `True` to run the database search
                in scan mode.

        Returns:
            `~pyncbitk.algo.SearchResultsSet`: The list of search results,
            with one `~pyncbitk.algo.SearchResults` item per query.

        """
        cdef TSeqLocVector         _queries_loc
        cdef TSeqLocVector         _subjects_loc
        cdef BlastSeqLoc           seqloc
        cdef CRef[IQueryFactory]   query_factory
        cdef CRef[IQueryFactory]   subject_factory
        cdef CRef[CLocalDbAdapter] db
        cdef CRef[CLocalBlast]     blast

        # prepare queries: a list of `BlastSeqLoc` objects
        if BlastQueries is BioSeq:
            if queries._ref.GetNonNullPointer().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = queries.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            query_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq](queries._ref)))
        elif BlastQueries is BioSeqSet:
            query_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq_set](queries._ref)))
        else:
            if not is_iterable(queries):
                queries = (queries, )
            for seqloc in queries:
                _queries_loc.push_back(seqloc._seqloc)
            query_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_queries_loc))

        # prepare subjects: either a list of `BlastSeqLoc` objects, or a `DatabaseReader`
        if BlastSubjects is DatabaseReader:
            _ty = subjects._ref.GetNonNullPointer().GetSequenceType()
            if _ty == ESeqType.eProtein:
                search_database = new CSearchDatabase(string(), EMoleculeType.eBlastDbIsProtein)
            elif _ty == ESeqType.eNucleotide:
                search_database = new CSearchDatabase(string(), EMoleculeType.eBlastDbIsNucleotide)
            else:
                raise ValueError(f"invalid sequence type: {_ty!r}")
            search_database.SetSeqDb(subjects._ref)
            db.Reset(new CLocalDbAdapter(search_database[0]))
        elif BlastSubjects is BioSeq:
            if subjects._ref.GetNonNullPointer().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = subjects.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq](subjects._ref)))
            db.Reset(new CLocalDbAdapter(subject_factory, self._opt, scan_mode))
        elif BlastSubjects is BioSeqSet:
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq_set](subjects._ref)))
            db.Reset(new CLocalDbAdapter(subject_factory, self._opt, scan_mode))
        else:
            if not is_iterable(subjects):
                subjects = (subjects, )
            for seqloc in subjects:
                _subjects_loc.push_back(seqloc._seqloc)
            subject_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_subjects_loc))
            db.Reset(new CLocalDbAdapter(subject_factory, self._opt, scan_mode))

        # prepare the BLAST program
        try:
            blast.Reset(new CLocalBlast(query_factory, self._opt, db))
        except RuntimeError as err:
            raise ValueError("Failed initializing BLAST") from err
        # if (m_InterruptFnx != NULL) {
        #     m_Blast->SetInterruptCallback(m_InterruptFnx, m_InterruptUserData);
        # }
        # // Set the hitlist size to the total number of subject sequences, to
        # // make sure that no hits are discarded (ported from CBl2Seq::SetupSearch
        # m_OptsHandle.SetHitlistSize((int) m_tSubjects.size());

        # run BLAST and get results
        with nogil:
            results = blast.GetNonNullPointer().Run()

        # check for warnings or errors
        messages = blast.GetNonNullPointer().GetSearchMessages()
        if messages.HasMessages():
            print(messages.ToString().decode())

        return SearchResultsSet._wrap(results)
        # messages = blast.GetSearchMessages() # TODO


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
        **kwargs,
    ):
        cdef CBlastAdvancedProteinOptionsHandle* handle = new CBlastAdvancedProteinOptionsHandle()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(**kwargs)


cdef class BlastN(NucleotideBlast):
    """A command object for running ``blastn`` searches.
    """

    def __init__(
        self,
        *,
        **kwargs,
    ):
        cdef CBlastNucleotideOptionsHandle* handle = new CBlastNucleotideOptionsHandle()
        handle.SetTraditionalBlastnDefaults()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(**kwargs)

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