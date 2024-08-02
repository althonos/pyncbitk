# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.utility cimport move

from .toolkit.algo.blast.api.bl2seq cimport CBl2Seq
from .toolkit.algo.blast.api.blast_types cimport EProgram, ProgramNameToEnum, TSeqAlignVector, EProgramToTaskName
from .toolkit.algo.blast.api.sseqloc cimport SSeqLoc, TSeqLocVector
from .toolkit.algo.blast.api.local_blast cimport CLocalBlast
from .toolkit.algo.blast.api.blast_options_handle cimport CBlastOptionsHandle, CBlastOptionsFactory
from .toolkit.algo.blast.api.blast_options_handle cimport CBlastOptionsHandle, CBlastOptionsFactory
from .toolkit.algo.blast.api.blast_nucl_options cimport CBlastNucleotideOptionsHandle
from .toolkit.algo.blast.api.query_data cimport IQueryFactory
from .toolkit.algo.blast.api.objmgr_query_data cimport CObjMgr_QueryFactory
from .toolkit.algo.blast.api.objmgrfree_query_data cimport CObjMgrFree_QueryFactory
from .toolkit.algo.blast.api.local_db_adapter cimport CLocalDbAdapter
from .toolkit.algo.blast.api.blast_results cimport CSearchResultSet, CSearchResults, size_type as CSearchResults_size_type
from .toolkit.corelib.ncbiobj cimport CConstRef, CRef
from .toolkit.corelib.ncbistr cimport kEmptyStr
from .toolkit.corelib.ncbitype cimport Uint4
from .toolkit.corelib.tempstr cimport CTempString
from .toolkit.objects.blastdb.blast_def_line_set cimport CBlast_def_line_set
from .toolkit.objects.general.object_id cimport CObject_id
from .toolkit.objects.general.object_id cimport E_Choice as CObject_id_choice
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objects.seq.bioseq cimport ELabelType as CBioseq_ELabelType
from .toolkit.objects.seq.seq_data cimport CSeq_data
from .toolkit.objects.seq.seq_data cimport E_Choice as CSeq_data_choice
from .toolkit.objects.seq.seq_inst cimport CSeq_inst, ETopology, EStrand
from .toolkit.objects.seq.seq_inst cimport EMol as CSeq_inst_mol
from .toolkit.objects.seq.seq_inst cimport ERepr as CSeq_inst_repr
from .toolkit.objects.seq.iupacna cimport CIUPACna
from .toolkit.objects.seqloc.seq_id cimport CSeq_id
from .toolkit.objects.seqloc.seq_id cimport E_Choice as CSeq_id_choice
from .toolkit.objects.seqloc.seq_loc cimport CSeq_loc
from .toolkit.objects.seqloc.seq_loc cimport E_Choice as CSeq_loc_choice
from .toolkit.objects.seqloc.textseq_id cimport CTextseq_id
from .toolkit.objects.seqset.seq_entry cimport CSeq_entry
from .toolkit.objects.seqset.seq_entry cimport E_Choice as CSeq_entry_choice
from .toolkit.objects.seqalign.seq_align cimport CSeq_align
from .toolkit.objects.seqalign.seq_align_set cimport CSeq_align_set
from .toolkit.objects.seqalign.score cimport CScore, C_Value as CScore_value, E_Choice as CScore_value_choice
from .toolkit.objmgr.object_manager cimport CObjectManager
from .toolkit.objmgr.scope cimport CScope
from .toolkit.objtools.readers.fasta cimport CFastaReader
from .toolkit.serial.serialbase cimport CSerialObject, MSerial_Format_AsnText
from .toolkit.serial.serialdef cimport ESerialRecursionMode
from .toolkit.objects.blastdb.blast_def_line cimport CBlast_def_line
from .toolkit.objects.blastdb.blast_def_line_set cimport CBlast_def_line_set
from .toolkit.algo.blast.api.uniform_search cimport CSearchDatabase, EMoleculeType
from .toolkit.objtools.blast.seqdb_reader.seqdb cimport CSeqDB, ESeqType

from .objects cimport ObjectId, SeqLoc, SeqAlignSet, SeqAlign, BioSeq
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

# --- BLAST results ------------------------------------------------------------

cdef class SearchResultsSet:
    cdef CRef[CSearchResultSet] _ref

    @staticmethod
    cdef SearchResultsSet _wrap(CRef[CSearchResultSet] ref):
        cdef SearchResultsSet obj = SearchResultsSet.__new__(SearchResultsSet)
        obj._ref = ref
        return obj

    def __len__(self):
        return self._ref.GetObject().size()

    def __getitem__(self, ssize_t index):
        cdef ssize_t _index  = index
        cdef ssize_t _length = self._ref.GetObject().size()

        if _index < 0:
            _index += _length
        if _index < 0 or _index >= _length:
            raise IndexError(index)

        cdef CSearchResultSet* obj  = &self._ref.GetObject()
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
    def alignments(self):
        return SeqAlignSet._wrap(self._ref.GetObject().GetSeqAlignMut())


# --- BLAST --------------------------------------------------------------------

cdef class Blast:
    cdef CRef[CBlastOptionsHandle] _opt

    @staticmethod
    def tasks():
        return [ x.decode() for x in CBlastOptionsFactory.GetTasks() ]

    def __init__(
        self,
        *,
        int window_size = 40,
    ):
        if self._opt.Empty():
            raise TypeError("Cannot instantiate abstract class Blast")
        self.window_size = window_size

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.program!r}, window_size={self.window_size})"

    # --- Properties -----------------------------------------------------------

    @property
    def program(self):
        cdef EProgram program = self._opt.GetObject().GetOptions().GetProgram()
        return EProgramToTaskName(program).decode('ascii')

    @property
    def window_size(self):
        """`int`: The window size for multiple hits.
        """
        return self._opt.GetObject().GetWindowSize()

    @window_size.setter
    def window_size(self, int window_size):
        if window_size < 0:
            raise ValueError(f"window_size must be a positive integer, got {window_size!r}")
        self._opt.GetObject().SetWindowSize(window_size)

    @property
    def off_diagonal_range(self):
        """`int`: The number of off-diagonals to search for the second hit.
        """
        return self._opt.GetObject().GetOffDiagonalRange()

    @property
    def xdrop_gap(self):
        """`float`: X-dropoff value (in bits) for preliminary gapped extensions.
        """
        return self._opt.GetObject().GetGapXDropoff()

    @property
    def xdrop_gap_final(self):
        """`float`: X-dropoff value (in bits) for final gapped alignment.
        """
        return self._opt.GetObject().GetGapXDropoffFinal()

    # @property
    # def max_hsps(self):
    #     """`int`: Maximum number of HSPs per subject to save for each query.
    #     """
    #     return self._opt.GetObject().GetMaxHspsPerSubject()

    @property
    def evalue(self):
        """`float`: Expectation value (E) threshold for saving hits.
        """
        return self._opt.GetObject().GetEvalueThreshold()

    @property
    def percent_identity(self):
        """`float`: Percentage identity threshold for saving hits.
        """
        return self._opt.GetObject().GetPercentIdentity()

    @property
    def coverage_hsp(self):
        """`float`: Query coverage percentage per HSP.
        """
        return self._opt.GetObject().GetQueryCovHspPerc()

    @property
    def gapped(self):
        """`bool`: `False` if alignments are performed in ungapped mode only.
        """
        return self._opt.GetObject().GetGappedMode()

    @property
    def culling_limit(self):
        """`int`: The culling limit for hits.

        If the query range of a hit is enveloped by that of at least this many
        higher-scoring hits, delete the hit.

        """
        return self._opt.GetObject().GetCullingLimit()

    @property
    def db_size(self):
        """`int`: The effective length of the database.
        """
        return self._opt.GetObject().GetDbLength()

    @property
    def search_space(self):
        """`int`: The effective length of the search space.
        """
        return self._opt.GetObject().GetEffectiveSearchSpace()

    # --- Public Methods -------------------------------------------------------

    cpdef SearchResultsSet run(
        self,
        queries,
        subjects,
        bool scan_mode = False
    ):
        cdef TSeqLocVector         _queries_loc
        cdef TSeqLocVector         _subjects_loc
        cdef BioSeq                _query_seq
        cdef BioSeq                _subject_seq
        cdef DatabaseReader        _db
        cdef BlastSeqLoc           seqloc
        cdef CRef[IQueryFactory]   query_factory
        cdef CRef[IQueryFactory]   subject_factory
        cdef CRef[CLocalDbAdapter] db
        cdef CRef[CLocalBlast]     blast

        # prepare queries: a list of `BlastSeqLoc` objects 
        if isinstance(queries, BioSeq):
            _query_seq = queries
            if _query_seq._ref.GetObject().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = _query_seq.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            query_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq](_query_seq._ref)))
        else:
            if not is_iterable(queries):
                queries = (queries, )
            for seqloc in queries:
                _queries_loc.push_back(seqloc._seqloc)
            query_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_queries_loc))

        # prepare subjects: either a list of `BlastSeqLoc` objects, or a `DatabaseReader`
        if isinstance(subjects, DatabaseReader):
            _db = subjects
            _ty = _db._ref.GetObject().GetSequenceType()
            if _ty == ESeqType.eProtein:
                search_database = new CSearchDatabase(string(), EMoleculeType.eBlastDbIsProtein)
            elif _ty == ESeqType.eNucleotide:
                search_database = new CSearchDatabase(string(), EMoleculeType.eBlastDbIsNucleotide)
            else:
                raise ValueError(f"invalid sequence type: {_ty!r}")
            search_database.SetSeqDb((<DatabaseReader> subjects)._ref)
            db.Reset(new CLocalDbAdapter(search_database[0]))
        elif isinstance(subjects, BioSeq):
            _subject_seq = subjects
            if _subject_seq._ref.GetObject().GetInst().GetRepr() != CSeq_inst_repr.eRepr_raw:
                ty = _subject_seq.instance.__class__.__name__
                raise ValueError(f"Unsupported instance type: {ty}")
            subject_factory.Reset(<IQueryFactory*> new CObjMgrFree_QueryFactory(CConstRef[CBioseq](_subject_seq._ref)))
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
            results = blast.GetObject().Run()
        return SearchResultsSet._wrap(results)
        # messages = blast.GetSearchMessages() # TODO


cdef class NucleotideBlast(Blast):
    pass

cdef class ProteinBlast(Blast):
    pass

cdef class MappingBlast(Blast):
    pass


cdef class BlastN(NucleotideBlast):

    def __init__(
        self,
        *,
        int window_size = 40,
    ):
        cdef CBlastNucleotideOptionsHandle* handle = new CBlastNucleotideOptionsHandle()
        handle.SetTraditionalBlastnDefaults()
        self._opt.Reset(<CBlastOptionsHandle*> handle)
        super().__init__(
            window_size=window_size,
        )
