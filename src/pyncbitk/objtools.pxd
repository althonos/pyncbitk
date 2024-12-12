# cython: language_level=3, linetrace=True, binding=True

from cython.operator cimport preincrement

from libc.limits cimport UINT_MAX
from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string
from libcpp.utility cimport move

from iostream cimport istream, streambuf

from .toolkit.corelib.ncbiobj cimport CRef
from .toolkit.corelib.ncbimisc cimport TSeqPos
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objtools.readers.fasta cimport CFastaReader
from .toolkit.objtools.blast.seqdb_reader.seqdb cimport CSeqDB, CSeqDBIter
from .toolkit.objtools.blast.seqdb_writer.writedb cimport CWriteDB
from .toolkit.objtools.alnmgr.alnmap cimport CAlnMap

from .objects.seqset cimport Entry
from .objects.seq cimport BioSeq
from .objects.seqalign cimport DenseSegments

import os

# --- FastaReader --------------------------------------------------------------

cdef class FastaReader:
    cdef CFastaReader* _reader
    cdef streambuf*    _buffer
    cdef istream*      _stream

    cpdef BioSeq read(self)

# --- BlastDatabase ------------------------------------------------------------

cdef class _DatabaseIter:
    cdef CRef[CSeqDB]   _ref
    cdef DatabaseReader db
    cdef CSeqDBIter*    it
    cdef size_t         length

cdef class DatabaseKeysIter(_DatabaseIter):
    pass

cdef class DatabaseValuesIter(_DatabaseIter):
    pass

cdef class DatabaseKeys:
    cdef DatabaseReader db
    cdef CRef[CSeqDB]   _ref

cdef class DatabaseValues:
    cdef DatabaseReader db
    cdef CRef[CSeqDB]   _ref

cdef class DatabaseReader:
    cdef CRef[CSeqDB] _ref

    cpdef DatabaseKeys keys(self)
    cpdef DatabaseValues values(self)

cdef class DatabaseWriter:
    cdef          CRef[CWriteDB] _ref
    cdef readonly bool           closed

# --- AlignMap -----------------------------------------------------------------

cdef class AlignMapRow:
    cdef          CRef[CAlnMap] _ref
    cdef          size_t        _index
    cdef readonly AlignMap      map

cdef class AlignMap:
    cdef          CRef[CAlnMap] _ref
    cdef readonly DenseSegments segments
