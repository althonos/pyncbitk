from libcpp cimport bool
from libcpp.string cimport string

from ...corelib.tempstr cimport CTempString
from ...corelib.ncbistr cimport kEmptyStr
from ...serial.serialbase cimport CSerialObject
from ...serial.serializable cimport CSerializable

cdef extern from "objects/seqloc/Textseq_id_.hpp" namespace "ncbi::objects::CTextseq_id_Base" nogil:
    ctypedef string TName
    ctypedef string TAccession
    ctypedef string TRelease
    ctypedef int TVersion
    
cdef extern from "objects/seqloc/Textseq_id_.hpp" namespace "ncbi::objects" nogil:

    cppclass CTextseq_id_Base(CSerialObject):

        CTextseq_id_Base()

        # optional
        bool IsSetName() const
        bool CanGetName() const
        void ResetName()
        const TName& GetName() const
        void SetName(const TName& value)
        void SetName(TName& value)
        TName& GetNameMut "SetName"()

        # optional
        bool IsSetAccession() const
        bool CanGetAccession() const
        void ResetAccession()
        const TAccession& GetAccession() const
        void SetAccession(const TAccession& value)
        void SetAccession(TAccession& value)
        TAccession& GetAccessionMut "SetAccession"()

        # optional
        bool IsSetRelease() const
        bool CanGetRelease() const
        void ResetRelease()
        const TRelease& GetRelease() const
        void SetRelease(const TRelease& value)
        void SetRelease(TRelease& value)
        TRelease& GetReleaseMut "SetRelease" ()

        # optional
        bool IsSetVersion() const
        bool CanGetVersion() const
        void ResetVersion()
        TVersion GetVersion() const
        void SetVersion(TVersion value)
        TVersion& GetVersionMut "SetVersion" ()


        void Reset()


cdef extern from "objects/seqloc/Textseq_id.hpp" namespace "ncbi::objects" nogil:

    cppclass CTextseq_id(CTextseq_id_Base, CSerializable):
        CTextseq_id()

        CTextseq_id& Set(
            const CTempString& acc_in,
            const CTempString& name_id           = kEmptyStr,
            int                version           = 0,
            const CTempString& release_in        = kEmptyStr,
            bool               allow_dot_version = true
        )