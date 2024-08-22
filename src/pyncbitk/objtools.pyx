# cython: language_level=3, linetrace=True, binding=True

from cython.operator cimport preincrement

from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport move

from .toolkit.corelib.ncbiobj cimport CRef
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objects.seqset.seq_entry cimport CSeq_entry
from .toolkit.objtools.readers.fasta cimport CFastaReader, EFlags as CFastaReader_Flags
from .toolkit.objtools.blast.seqdb_reader.seqdb cimport CSeqDB, CSeqDBIter, ESeqType
from .toolkit.objtools.blast.seqdb_reader.seqdbcommon cimport EBlastDbVersion, EOidMaskType
from .toolkit.objtools.blast.seqdb_writer.writedb cimport CWriteDB, EIndexType

from .objects.seqset cimport Entry
from .objects.seq cimport BioSeq

import os

# --- FastaReader --------------------------------------------------------------

cdef class FastaReader:
    """An iterative reader over the sequences in a FASTA file.
    """

    def __cinit__(self):
        self._reader = NULL

    def __dealloc__(self):
        del self._reader

    def __init__(
        self,
        object path,
        *,
        bool split = True,
    ):
        cdef bytes _path = os.fsencode(path)
        cdef int   flags = 0

        if not split:
            flags |= CFastaReader_Flags.fNoSplit

        self._reader = new CFastaReader(_path, flags)

    def __iter__(self):
        return self

    def __next__(self):
        cdef BioSeq seq = self.read()
        if seq is None:
            raise StopIteration
        return seq

    cpdef BioSeq read(self):
        assert self._reader != NULL

        cdef CRef[CSeq_entry] entry
        cdef CBioseq*         bioseq

        if self._reader.AtEOF():
            return None

        entry = self._reader.ReadOneSeq()
        bioseq = &entry.GetNonNullPointer().GetSeqMut()
        return BioSeq._wrap(CRef[CBioseq](bioseq))

# --- BlastDatabase ------------------------------------------------------------

cdef class DatabaseIter:
    """An iterator over the sequences of a BLAST database.
    """

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
    """A handle allowing to read the contents of a BLAST database.
    """

    @staticmethod
    def _search_path():
        return CSeqDB.GenerateSearchPath().decode()

    def __init__(self, object name, str type = None):
        """__init__(self, name, type=None)\n--\n

        Create a new reader for a database of the given name.

        Arguments:
            name (`str`): The name of the database, as given when the
                database was created.
            type (`str` or `None`): The type of sequences in the database.
                If `None` given, the database type will be detected from
                the metadata.

        """
        cdef bytes   _name = name.encode()  # FIXME: os.fsencode?
        cdef CSeqDB* _db   = new CSeqDB(<string> _name, ESeqType.eUnknown)
        self._ref.Reset(_db)

    def __iter__(self):
        return DatabaseIter(self)

    def __len__(self):
        return self._ref.GetObject().GetNumSeqs()


cdef class DatabaseWriter:
    """A handle allowing to write sequences to a BLAST database.
    """

    def __init__(
        self,
        name,
        str type not None = "nucleotide",
        *,
        object title = None,
        int version = 4,

    ):
        """__init__(self, name, type="nucleotide", *, title=None, version=4)\n--\n

        Create a new database writer.

        Arguments:
            name (`str`): The name of the database, which is used as a path
                prefix to create the database files.
            type (`str`): Either ``nucleotide`` for a nucleotide database,
                or ``"protein"`` for a protein database.
            title (`str` or `None`): The title of the database.
            version (`int`): The database format version, either ``4`` (the 
                default) or ``5``.

        """
        cdef bytes           _path
        cdef bytes           _title
        cdef ESeqType        dbtype
        cdef EBlastDbVersion dbver
        cdef CWriteDB*       writer

        _path = os.fsencode(name)
        _parent = os.path.dirname(_path)
        if _parent and not os.path.exists(_parent):
            raise FileNotFoundError(os.fsdecode(_parent))

        if type == "nucleotide":
            dbtype = ESeqType.eNucleotide
        elif type == "protein":
            dbtype = ESeqType.eProtein
        else:
            raise ValueError(f"type must be either 'nucleotide' or 'protein', got {type!r}")

        if version == 4:
            dbver = EBlastDbVersion.eBDB_Version4
        elif version == 5:
            dbver = EBlastDbVersion.eBDB_Version5
        else:
            raise ValueError(f"version must be either 4 or 5, got {version!r}")

        if title is None:
            _title = _path
        elif isinstance(title, str):
            _title = title.encode()
        else:
            _title = title

        writer = new CWriteDB(
            _path,
            dbtype,
            _title,
            EIndexType.eDefault,
            True,
            False,
            False,
            dbver,
            False,
            EOidMaskType.fNone,
            False,
        )

        self._ref.Reset(writer)
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return None

    @property
    def volumes(self):
        """`list` or `str`: The list of volumes written by the writer.
        """
        cdef vector[string] volumes
        self._ref.GetNonNullPointer().ListVolumes(volumes)
        return [os.fsdecode(v) for v in volumes]

    @property
    def files(self):
        """`list` or `str`: The list of files written by the writer.
        """
        cdef vector[string] files
        self._ref.GetNonNullPointer().ListFiles(files)
        return [os.fsdecode(f) for f in files]

    def append(self, BioSeq sequence):
        """Add a sequence to the database.

        Arguments:
            sequence (`~pyncbitk.objects.seq.BioSeq`): The sequence to add
                to the database.

        """
        cdef CWriteDB* writer = self._ref.GetNonNullPointer()
        writer.AddSequence(sequence._ref.GetObject())

    def close(self):
        """Close the database and write the remaining buffered sequences.
        """
        if not self.closed:
            self._ref.GetNonNullPointer().Close()
            self.closed = True

    