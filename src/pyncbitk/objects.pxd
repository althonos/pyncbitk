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
from .toolkit.objects.seqset.bioseq_set cimport CBioseq_set
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

# --- ObjectId -----------------------------------------------------------------

cdef class ObjectId:
    cdef CRef[CObject_id] _ref

    @staticmethod
    cdef ObjectId _wrap(CRef[CObject_id] ref)
     
cdef class StrId(ObjectId):
    pass

cdef class IntId(ObjectId):
    pass

# --- TextSeqId ----------------------------------------------------------------

cdef class TextSeqId:
    cdef CRef[CTextseq_id] _ref

    @staticmethod
    cdef TextSeqId _wrap(CRef[CTextseq_id] ref)

# --- SeqId --------------------------------------------------------------------

cdef class SeqId:
    cdef CRef[CSeq_id] _ref

    @staticmethod
    cdef SeqId _wrap(CRef[CSeq_id] ref)


cdef class LocalId(SeqId):
    pass
        
cdef class RefSeqId(SeqId):
    pass

cdef class GenBankId(SeqId):
    pass
    
cdef class ProteinDataBankId(SeqId):
    pass

cdef class GeneralId(SeqId):
    pass

# --- BioSeq -------------------------------------------------------------------

cdef class BioSeq:
    cdef CRef[CBioseq] _ref

    @staticmethod
    cdef BioSeq _wrap(CRef[CBioseq] ref)
        

cdef class BioSeqSet:
    cdef CRef[CBioseq_set] _ref


# --- SeqInst ------------------------------------------------------------------

cdef class SeqInst:
    cdef CRef[CSeq_inst] _ref

    @staticmethod
    cdef SeqInst _wrap(CRef[CSeq_inst] ref)

cdef class EmptyInst(SeqInst):
    pass

cdef class VirtualInst(SeqInst):
    pass

cdef class ContinuousInst(SeqInst):
    pass

cdef class SegmentedInst(SeqInst):
    pass

cdef class ConstructedInst(SeqInst):
    pass

cdef class RefInst(SeqInst):
    pass

cdef class ConsensusInst(SeqInst):
    pass

cdef class MapInst(SeqInst):
    pass

cdef class DeltaInst(SeqInst):
    
    cpdef ContinuousInst to_continuous(self)


# --- SeqData ------------------------------------------------------------------

cdef class SeqData:
    cdef CRef[CSeq_data] _ref

    @staticmethod
    cdef SeqData _wrap(CRef[CSeq_data] ref)

cdef class SeqAaData(SeqData):
    pass

cdef class SeqNaData(SeqData):
    pass

cdef class IupacNaData(SeqNaData):
    cpdef str decode(self)


cdef class IupacAaData(SeqAaData):
    pass

cdef class Ncbi2NaData(SeqNaData):
    pass

cdef class Ncbi4NaData(SeqNaData):
    pass

cdef class Ncbi8NaData(SeqNaData):
    pass

cdef class NcbiPNaData(SeqNaData):
    pass

cdef class Ncbi8AaData(SeqAaData):
    pass

cdef class NcbiEAaData(SeqAaData):
    pass

cdef class NcbiPAaData(SeqAaData):
    pass

cdef class NcbiStdAa(SeqAaData):
    pass

cdef class GapData(SeqData):
    pass


# --- Entry --------------------------------------------------------------------

cdef class Entry:
    cdef CRef[CSeq_entry] _ref

    @staticmethod
    cdef Entry _wrap(CRef[CSeq_entry] ref)

cdef class SeqEntry(Entry):
    pass

cdef class SetEntry(Entry):
    pass


# --- SeqLoc -------------------------------------------------------------------

cdef class SeqLoc:
    cdef CRef[CSeq_loc] _loc

cdef class NullLoc(SeqLoc):
    pass

cdef class EmptySeqLoc(SeqLoc):
    pass

cdef class WholeSeqLoc(SeqLoc):
    pass

cdef class SeqIntervalLoc(SeqLoc):
    pass

cdef class PackedSeqLoc(SeqLoc):
    pass

cdef class PackedPointsLoc(SeqLoc):
    pass

cdef class MixLoc(SeqLoc):
    pass

cdef class EquivalentLoc(SeqLoc):
    pass

cdef class BondLoc(SeqLoc):
    pass

cdef class FeatureLoc(SeqLoc):
    pass

# --- SeqAlign -----------------------------------------------------------------

cdef class SeqAlignScore:
    cdef CRef[CScore] _ref

    @staticmethod
    cdef SeqAlignScore _wrap(CRef[CScore] ref)

cdef class SeqAlign:
    cdef CRef[CSeq_align] _ref

    @staticmethod
    cdef SeqAlign _wrap(CRef[CSeq_align] ref)

cdef class SeqAlignSet:
    cdef CRef[CSeq_align_set] _ref

    @staticmethod
    cdef SeqAlignSet _wrap(CRef[CSeq_align_set] ref)