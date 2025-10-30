# cython: language_level=3, linetrace=True, binding=True

from cython.operator cimport preincrement

from libc.limits cimport UINT_MAX
from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport move

from iostream cimport istream

from .toolkit.corelib.ncbiobj cimport CRef
from .toolkit.corelib.ncbimisc cimport TSeqPos
from .toolkit.corelib.ncbistre cimport CNcbiIstream
from .toolkit.objects.seq.bioseq cimport CBioseq
from .toolkit.objects.seqset.seq_entry cimport CSeq_entry
from .toolkit.objtools.readers.fasta cimport CFastaReader, EFlags as CFastaReader_Flags
from .toolkit.objtools.blast.seqdb_reader.seqdb cimport CSeqDB, CSeqDBIter, ESeqType
from .toolkit.objtools.blast.seqdb_reader.seqdbcommon cimport EBlastDbVersion, EOidMaskType
from .toolkit.objtools.blast.seqdb_writer.writedb cimport CWriteDB, EIndexType
from .toolkit.objtools.alnmgr.alnmap cimport CAlnMap, TNumrow

from pystreambuf cimport pyreadbuf, pyreadintobuf

from .objects.seqset cimport Entry
from .objects.seqid cimport SeqId
from .objects.seq cimport BioSeq
from .objects.seqalign cimport DenseSegments

import os
import sys

# --- FastaReader --------------------------------------------------------------

cdef class FastaReader:
    """An iterative reader over the sequences in a FASTA file.
    """

    def __cinit__(self):
        self._reader = NULL
        self._buffer = NULL
        self._stream = NULL

    def __dealloc__(self):
        del self._reader
        if self._stream is not NULL:
            del self._stream
            self._stream = NULL
        if self._buffer is not NULL:
            del self._buffer
            self._buffer = NULL

    def __init__(
        self,
        object file not None,
        *,
        bool split = True,
        bool parse_ids = True,
    ):
        """__init__(self, file, *, split=True, parse_ids=True)\n--\n

        Create a new FASTA reader from a file or a file-like object.

        Arguments:
            file (`str`, `os.PathLike` or file-like object): Either the path
                to a file to be open, or a Python file-like object open in
                binary mode.
            split (`bool`): Set to `False` to force the reader to produce
                `~pyncbitk.objects.seq.BioSeq` objects where the instance
                is a `~pyncbitk.objects.seqinst.ContinuousInst` object.
            parse_ids (`bool`): If `True` (the default), the FASTA header
                line will be parsed to build a sequence identifier. This
                may cause errors on certain files where the IDs are too
                long (>50 letters). Pass `False` to generate a `LocalId`
                for each `BioSeq` using a counter. The whole FASTA header
                line can be accessed as `~pyncbitk.objects.seqdesc.TitleDesc`
                instances of the `BioSeq` in both cases.

        """
        cdef bytes      path
        cdef int        flags = 0

        if not split:
            flags |= CFastaReader_Flags.fNoSplit
        if not parse_ids:
            flags |= CFastaReader_Flags.fNoParseID

        if hasattr(file, "readinto") and sys.implementation.name == "cpython":
            self._buffer = new pyreadintobuf(file)
            self._stream = new istream(self._buffer)
            self._reader = new CFastaReader(self._stream[0], flags)
        elif hasattr(file, "read"):
            self._buffer = new pyreadbuf(file)
            self._stream = new istream(self._buffer)
            self._reader = new CFastaReader(self._stream[0], flags)
        else:
            path = os.fsencode(file)
            self._reader = new CFastaReader(path, flags)

    def __iter__(self):
        return self

    def __next__(self):
        cdef BioSeq seq = self.read()
        if seq is None:
            raise StopIteration
        return seq

    @property
    def max_id_length(self):
        return self._reader.GetMaxIDLength()

    cpdef BioSeq read(self):
        """Read a single sequence if available.

        Returns:
            `~pyncbitk.objects.seq.BioSeq`: The next sequence in the FASTA
            file, or `None` if the reader reached the end of the file.

        """
        assert self._reader != NULL

        cdef CRef[CSeq_entry] entry
        cdef CBioseq*         bioseq

        if self._reader.AtEOF():
            return None

        entry = self._reader.ReadOneSeq()
        bioseq = &entry.GetNonNullPointer().GetSeqMut()
        return BioSeq._wrap(CRef[CBioseq](bioseq))

# --- BlastDatabase ------------------------------------------------------------

cdef class _DatabaseIter:

    def __init__(self, DatabaseReader db):
        if self.it is not NULL:
            del self.it
            self.it = NULL
        self.db = db
        self._ref = db._ref
        self.it = new CSeqDBIter(db._ref.GetNonNullPointer().Begin())
        self.length = self._ref.GetNonNullPointer().GetNumSeqs()

    def __dealloc__(self):
        del self.it
        self.it = NULL

    def __iter__(self):
        return self

    def __len__(self):
        return self.length

    def __next__(self):
        raise StopIteration

cdef class DatabaseKeysIter(_DatabaseIter):
    """An interator over the sequence identifiers of a BLAST database.
    """

    def __next__(self):
        cdef int           oid
        cdef CRef[CBioseq] seq

        if not self.it[0]:
            raise StopIteration

        oid = self.it.GetOID()
        ids = self._ref.GetNonNullPointer().GetSeqIDs(oid)
        preincrement(self.it[0])
        self.length -= 1
        return SeqId._wrap(ids.front())

cdef class DatabaseValuesIter(_DatabaseIter):
    """An iterator over the sequences of a BLAST database.
    """

    def __next__(self):
        cdef int           oid
        cdef CRef[CBioseq] seq

        if not self.it[0]:
            raise StopIteration

        oid = self.it.GetOID()
        seq = self._ref.GetNonNullPointer().GetBioseq(oid)
        preincrement(self.it[0])
        self.length -= 1
        return BioSeq._wrap(seq)

cdef class DatabaseKeys:
    """A set-like view over the keys of a BLAST database.
    """

    def __init__(self, DatabaseReader db not None):
        self.db = db
        self._ref = db._ref

    def __iter__(self):
        return DatabaseKeysIter(self.db)

    def __len__(self):
        return len(self.db)

    def __contains__(self, object item):
        cdef SeqId sid
        cdef int   oid
        if not isinstance(item, SeqId):
            return False
        sid = item
        return self._ref.GetObject().SeqidToOid(sid._ref.GetObject(), oid)


cdef class DatabaseValues:
    """A set-like view over the values of a BLAST database.
    """

    def __init__(self, DatabaseReader db not None):
        self.db = db
        self._ref = db._ref

    def __iter__(self):
        return DatabaseValuesIter(self.db)

    def __len__(self):
        return len(self.db)

    def __contains__(self, object item):
        cdef SeqId         sid
        cdef int           oid
        cdef CRef[CBioseq] seq
        if not isinstance(item, BioSeq):
            return False
        sid = item.id
        if not self._ref.GetObject().SeqidToOid(sid._ref.GetObject(), oid):
            return False
        seq = self._ref.GetNonNullPointer().GetBioseq(oid)
        return BioSeq._wrap(seq) == item


cdef class DatabaseReader:
    """A reader over the contents of a BLAST database.

    This class implements access to the BLAST database as a
    `~collections.abc.Mapping` of `~pyncbitk.objects.seqid.SeqId` objects
    to `~pyncbitk.objects.seq.BioSeq` objects.

    """

    @staticmethod
    def _search_path():
        return CSeqDB.GenerateSearchPath().decode()

    def __init__(
        self,
        object name not None,
        str type = None
    ):
        """__init__(self, name, type=None)\n--\n

        Create a new reader for a database of the given name.

        Arguments:
            name (`str` or `os.PathLike`): The name of the database, as given
                when the database was created.
            type (`str` or `None`): The type of sequences in the database,
                either ``nucleotide`` or ``protein``. If `None` given,
                the database type will be detected from the metadata.

        """
        # TODO: handle type given in argument
        cdef bytes   _name =  os.fsencode(name)
        cdef CSeqDB* _db   = new CSeqDB(<string> _name, ESeqType.eUnknown)
        self._ref.Reset(_db)

    def __iter__(self):
        return DatabaseKeysIter(self)

    def __len__(self):
        return self._ref.GetNonNullPointer().GetNumSeqs()

    def __getitem__(self, object index):
        cdef SeqId         seq_id
        cdef CRef[CBioseq] bioseq

        if not isinstance(index, SeqId):
            ty = type(index).__name__
            raise TypeError("database keys must be SeqId, not {ty!r}")

        seq_id = index
        bioseq = self._ref.GetNonNullPointer().SeqidToBioseq(seq_id._ref.GetObject())
        return BioSeq._wrap(bioseq)

    @property
    def version(self):
        """`int`: The database format version.
        """
        cdef EBlastDbVersion dbver = self._ref.GetNonNullPointer().GetBlastDbVersion()
        if dbver == EBlastDbVersion.eBDB_Version4:
            return 4
        elif dbver == EBlastDbVersion.eBDB_Version5:
            return 5
        raise RuntimeError()

    cpdef DatabaseKeys keys(self):
        """Get a set-like view over the keys of the database.

        Returns:
            `~pyncbitk.objtools.DatabaseKeys`: The keys of the database,
            i.e. the `~pyncbitk.objects.seqid.SeqId` of the sequences
            stored in the database.

        """
        return DatabaseKeys(self)

    cpdef DatabaseValues values(self):
        """Get a set-like view over the values of the database.

        Returns:
            `~pyncbitk.objtools.DatabaseValues`: The values of the database,
            i.e. the `~pyncbitk.objects.seq.BioSeq` storing the sequences
            in the database.

        """
        return DatabaseValues(self)


cdef class DatabaseWriter:
    """A handle allowing to write sequences to a BLAST database.
    """

    def __init__(
        self,
        object name not None,
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
            type (`str`): Either ``"nucleotide"`` for a nucleotide database,
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
        """`list` of `str`: The list of files written by the writer.
        """
        cdef vector[string] files
        self._ref.GetNonNullPointer().ListFiles(files)
        return [os.fsdecode(f) for f in files]

    def append(self, BioSeq sequence not None):
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


# --- AlignMap -----------------------------------------------------------------

cdef class AlignMap:
    """A helper class to handle coordinates of `DenseSegments` of an alignment.
    """

    def __init__(self, DenseSegments segments):
        self.segments = segments
        self._ref.Reset(new CAlnMap(segments._ref.GetObject().GetDenseg()))

    def __len__(self):
        return self._ref.GetObject().GetNumRows()

    def __getitem__(self, index):
        cdef AlignMapRow row
        cdef ssize_t     index_ = index
        cdef ssize_t     length = self._ref.GetObject().GetNumRows()

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)

        row = AlignMapRow.__new__(AlignMapRow)
        row.map = self
        row._ref = self._ref
        row._index = <size_t> index_
        return row


cdef class AlignMapRow:
    """A row and its coordinates in an `AlignMap`.
    """

    @property
    def align_start(self):
        return self._ref.GetObject().GetSeqAlnStart(self._index)

    @property
    def align_stop(self):
        return self._ref.GetObject().GetSeqAlnStop(self._index)

    @property
    def sequence_start(self):
        return self._ref.GetObject().GetSeqStart(self._index)

    @property
    def sequence_stop(self):
        return self._ref.GetObject().GetSeqStop(self._index)

    @property
    def strand(self):
        return self._ref.GetObject().StrandSign(self._index)