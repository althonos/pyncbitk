# cython: language_level=3

from libcpp.string cimport string

from ..toolkit.objects.taxon1.taxon1 cimport CTaxon1
from ..toolkit.objects.taxon1.taxon2_data cimport CTaxon2_data
from ..toolkit.serial.serialbase cimport CSerialObject
from ..toolkit.corelib.ncbiobj cimport CRef

from ..serial cimport Serial


cdef class TaxonomyClient:
    cdef CTaxon1 _client

    def __init__(
        self,
    ):
        if not self._client.Init():
            self._last_error()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    cdef object _last_error(self) except *:
        cdef string err = self._client.GetLastError()
        cdef bytes  e   = err
        if e.startswith(b"ERROR: "):
            e = e[7:]
        raise RuntimeError(e.decode())

    cpdef TaxonData get(self, int taxid):
        cdef CRef[CTaxon2_data] ref = self._client.GetById(taxid)
        if ref.IsNull():
            self._last_error()
        return TaxonData._wrap(ref)

    cpdef void close(self):
        self._client.Fini()        


cdef class TaxonData(Serial):
    cdef CRef[CTaxon2_data] _ref

    @staticmethod
    cdef TaxonData _wrap(CRef[CTaxon2_data] ref):
        cdef TaxonData data = TaxonData.__new__(TaxonData)
        data._ref = ref
        return data

    cdef CSerialObject* _serial(self):
        return <CSerialObject*> self._ref.GetNonNullPointer()