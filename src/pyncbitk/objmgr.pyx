# cython: language_level=3, linetrace=True, binding=True

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.vector cimport vector
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

from .objects cimport BioSeq


cdef class ObjectManager:
    def __init__(self):
        self._mgr = CObjectManager.GetInstance()

    cpdef Scope scope(self):
        return Scope(self)

    def register_data_loader(self, str name not None):
        cdef bytes _name = name.encode()
        self._mgr.GetObject().RegisterDataLoader(NULL, <string> _name)

    def get_registered_names(self):
        cdef vector[string] names
        self._mgr.GetObject().GetRegisteredNames(names)
        return list(names)


cdef class Scope:
    def __init__(self, ObjectManager manager):
        self._scope.Reset(new CScope(manager._mgr.GetObject()))

    def __enter__(self):
        return self

    def __exit__(self, exc_ty, exc_val, traceback):
        self.close()
        return False

    def close(self):
        self._scope.ReleaseOrNull()

    cpdef void add_bioseq(self, BioSeq seq) except *:
        if self._scope.Empty():
            raise RuntimeError("attempted to use a closed scope")

        self._scope.GetObject().AddBioseq(seq._ref.GetObject())