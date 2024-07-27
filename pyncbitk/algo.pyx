# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.utility cimport move

from .toolkit.algo.blast.api.bl2seq cimport CBl2Seq
from .toolkit.algo.blast.api.blast_types cimport EProgram, ProgramNameToEnum, TSeqAlignVector
from .toolkit.algo.blast.api.sseqloc cimport SSeqLoc, TSeqLocVector
from .toolkit.algo.blast.api.local_blast cimport CLocalBlast
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

import os

from .objects cimport ObjectId, SeqLoc
from .objmgr cimport Scope

# --- BLAST input --------------------------------------------------------------

cdef class BlastSeqLoc:
    cdef SSeqLoc _seqloc

    def __init__(self, SeqLoc loc, Scope scope, SeqLoc mask = None):
        self._seqloc = SSeqLoc(loc._loc.GetObject(), scope._scope.GetObject())


cdef class Blast:
    cdef CRef[CBl2Seq] _blast

    def __init__(
        self, 
        BlastSeqLoc query, 
        object subject, 
        str program
    ): 

        cdef EProgram      p
        cdef BlastSeqLoc   s
        cdef TSeqLocVector subjects = TSeqLocVector()
        cdef bytes         _program = program.encode()

        try:
            p = ProgramNameToEnum(_program)
        except Exception as e:
            raise ValueError(f"Invalid BLAST program: {program!r}")

        if isinstance(subject, BlastSeqLoc):
            subject = (subject, )

        for s in subject:
            subjects.push_back(s._seqloc)

        


        self._blast.Reset(
            new CBl2Seq(query._seqloc, subjects, p)
        )

    def run(self):
        cdef TSeqAlignVector _raw_alignments 
        cdef SeqAlignSet     ali
        cdef list            alignments      
        
        with nogil:
            _raw_alignments = self._blast.GetObject().Run()

        alignments = []
        for ref in _raw_alignments:
            alignments.append(SeqAlignSet._wrap(ref))

        return alignments


cdef class SeqAlignScore:
    # TODO: inherit
    cdef CRef[CScore] _ref

    @staticmethod
    cdef SeqAlignScore _wrap(CRef[CScore] ref):
        cdef SeqAlignScore score = SeqAlignScore.__new__(SeqAlignScore)
        score._ref = ref
        return score

    def __iter__(self):
        yield self.id
        yield self.value

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(id={self.id!r}, value={self.value!r})"

    @property
    def id(self):
        if not self._ref.GetObject().IsSetId():
            return None
        id_ = &self._ref.GetObject().GetIdRw()
        cref = CRef[CObject_id](id_)
        return ObjectId._wrap(cref)

    @property
    def value(self):
        value = &self._ref.GetObject().GetValueRw()
        kind = value.Which()
        if kind == CScore_value_choice.e_Int:
            return value.GetInt()
        elif kind == CScore_value_choice.e_Real:
            return value.GetReal()
        raise TypeError(f"Unknown value type: {kind!r}")


cdef class SeqAlign:
    cdef CRef[CSeq_align] _ref

    @staticmethod
    cdef SeqAlign _wrap(CRef[CSeq_align] ref):
        cdef SeqAlign obj = SeqAlign.__new__(SeqAlign)
        obj._ref = ref
        return obj

    def dumps(self):
        cdef string s = string()
        s << <CSerialObject&> self._ref.GetNonNullPointer()[0]
        return s.decode()

    @property
    def scores(self):
        cdef CRef[CScore]  ref
        cdef SeqAlignScore score 
        cdef list          scores = []

        if not self._ref.GetObject().IsSetScore():
            return None

        for ref in self._ref.GetObject().GetScoreRw():
            scores.append(SeqAlignScore._wrap(ref))

        return scores

cdef class SeqAlignSet:
    cdef CRef[CSeq_align_set] _ref

    @staticmethod
    cdef SeqAlignSet _wrap(CRef[CSeq_align_set] ref):
        cdef SeqAlignSet obj = SeqAlignSet.__new__(SeqAlignSet)
        obj._ref = ref
        return obj

    def __iter__(self):
        cdef CRef[CSeq_align] ref
        for ref in self._ref.GetObject().GetRw():
            yield SeqAlign._wrap(ref)

    def __len__(self):
        return self._ref.GetObject().Get().size()
