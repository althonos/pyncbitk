# cython: language_level=3, linetrace=True, binding=True

from cython.operator cimport preincrement

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.utility cimport move

from .toolkit.corelib.ncbiobj cimport CRef
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objtools.readers.fasta cimport CFastaReader
from .toolkit.objtools.blast.seqdb_reader.seqdb cimport CSeqDB, CSeqDBIter, ESeqType

from .objects cimport Entry, BioSeq

import os

# --- FastaReader --------------------------------------------------------------

cdef class FastaReader:
    cdef CFastaReader* _reader

    cpdef Entry read(self)

# --- BlastDatabase ------------------------------------------------------------

cdef class DatabaseIter:
    cdef CRef[CSeqDB]         _ref
    cdef DatabaseReader       db
    cdef CSeqDBIter*          it

cdef class DatabaseReader:
    cdef CRef[CSeqDB] _ref