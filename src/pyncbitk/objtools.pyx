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

    def __cinit__(self):
        self._reader = NULL

    def __dealloc__(self):
        del self._reader

    def __init__(self, path):
        cdef bytes _path = os.fsencode(path)
        self._reader = new CFastaReader(_path)

    def read(self):
        assert self._reader != NULL
        _entry = self._reader.ReadOneSeq()
        return Entry._wrap(_entry)

# --- BlastDatabase ------------------------------------------------------------

cdef class DatabaseIter:

    def __init__(self, DatabaseReader db):

        if self.it is not NULL:
            del self.it
            self.it = NULL

        self.db = db
        self._ref = db._ref
        self.it = new CSeqDBIter(db._ref.GetObject().Begin())

    def __dealloc__(self):
        del self.it
        self.it = NULL

    def __iter__(self):
        return self

    def __next__(self):
        cdef int           oid
        cdef CRef[CBioseq] seq

        if not self.it[0]:
            raise StopIteration

        oid = self.it.GetOID()
        seq = self._ref.GetObject().GetBioseq(oid)
        preincrement(self.it[0])
        return BioSeq._wrap(seq)


cdef class DatabaseReader:

    @staticmethod
    def _search_path():
        return CSeqDB.GenerateSearchPath().decode()

    def __init__(self, object name, str type = None):
        cdef bytes   _name = name.encode()  # FIXME: os.fsencode?
        cdef CSeqDB* _db   = new CSeqDB(<string> _name, ESeqType.eUnknown)
        self._ref.Reset(_db)

    def __iter__(self):
        return DatabaseIter(self)

    def __len__(self):
        return self._ref.GetObject().GetNumSeqs()