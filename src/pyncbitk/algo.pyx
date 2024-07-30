# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.utility cimport move

from .toolkit.algo.blast.api.bl2seq cimport CBl2Seq
from .toolkit.algo.blast.api.blast_types cimport EProgram, ProgramNameToEnum, TSeqAlignVector
from .toolkit.algo.blast.api.sseqloc cimport SSeqLoc, TSeqLocVector
from .toolkit.algo.blast.api.local_blast cimport CLocalBlast
from .toolkit.algo.blast.api.blast_options_handle cimport CBlastOptionsHandle, CBlastOptionsFactory
from .toolkit.algo.blast.api.query_data cimport IQueryFactory
from .toolkit.algo.blast.api.objmgr_query_data cimport CObjMgr_QueryFactory
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

from .objects cimport ObjectId, SeqLoc, SeqAlignSet, SeqAlign
from .objmgr cimport Scope

import os


cdef extern from * nogil:
    """
    template <typename T>
    std::string dump(ncbi::CSerialObject& source) {
        std::string s;
        s << source;
        return s;
    }
    """
    cdef string dump[T](T source)


# --- BLAST input --------------------------------------------------------------

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

    def __init__(self, str program):
        cdef EProgram p
        cdef CBlastOptionsHandle* opt 
        cdef bytes    _program = program.encode()

        try:
            p = ProgramNameToEnum(_program)
        except Exception as e:
            raise ValueError(f"Invalid BLAST program: {program!r}") from e

        opt = CBlastOptionsFactory.Create(p)
        self._opt.Reset(opt)

    def run(self, queries, subjects, bool scan_mode = False):
        
        cdef TSeqLocVector         _queries  
        cdef TSeqLocVector         _subjects
        cdef BlastSeqLoc           seqloc
        cdef CRef[IQueryFactory]   query_factory      
        cdef CRef[IQueryFactory]   subject_factory    
        cdef CRef[CLocalDbAdapter] db             
        cdef CRef[CLocalBlast]     blast

        if isinstance(queries, BlastSeqLoc):
            queries = (queries, )
        if isinstance(subjects, BlastSeqLoc):
            subjects = (subjects, )

        for seqloc in queries:
            _queries.push_back(seqloc._seqloc)
        for seqloc in subjects:
            _subjects.push_back(seqloc._seqloc)

        query_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_queries))
        subject_factory.Reset(<IQueryFactory*> new CObjMgr_QueryFactory(_subjects))
        db.Reset(new CLocalDbAdapter(subject_factory, self._opt, scan_mode))
        blast.Reset(new CLocalBlast(query_factory, self._opt, db))
        # if (m_InterruptFnx != NULL) {
        #     m_Blast->SetInterruptCallback(m_InterruptFnx, m_InterruptUserData);
        # }
        # // Set the hitlist size to the total number of subject sequences, to 
        # // make sure that no hits are discarded (ported from CBl2Seq::SetupSearch
        # m_OptsHandle.SetHitlistSize((int) m_tSubjects.size()); 

        results = blast.GetObject().Run()
        return SearchResultsSet._wrap(results)
        # messages = blast.GetSearchMessages() # TODO


