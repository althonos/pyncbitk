# cython: language_level=3, linetrace=True, binding=True

from libcpp.string cimport string

from .toolkit.corelib.ncbistre cimport ostream, stringbuf, CNcbiOstream
from .toolkit.serial.serialbase cimport CSerialObject, MSerial_Format, MSerial_Flags, TSerial_Format_Flags
from .toolkit.serial.serialdef cimport ESerialRecursionMode, ESerialDataFormat, ESerial_Xml_Flags


# NOTE: Cython is not happy with the way the object serialization mechanism
#       is supposed to work because of const references; making it a simple
#       inline function is enough here.
cdef extern from * nogil:
    """
    void DumpObject(ncbi::CNcbiOstream& os, ncbi::CSerialObject& obj, ncbi::MSerial_Format& fmt) {
        os << fmt << obj;
    }
    """
    void DumpObject(CNcbiOstream& os, CSerialObject& obj, MSerial_Format& fmt) except +

cdef dict _SERIAL_DATAFORMAT_STR = {
    ESerialDataFormat.eSerial_None: None,
    ESerialDataFormat.eSerial_AsnText: "asntext",
    ESerialDataFormat.eSerial_AsnBinary: "asnbinary",
    ESerialDataFormat.eSerial_Xml: "xml",
    ESerialDataFormat.eSerial_Json: "json",
}

cdef dict _SERIAL_DATAFORMAT_ENUM = {
    v:k for k,v in _SERIAL_DATAFORMAT_STR.items()
}

cdef class Serial:
    """Abstract base class for objects part of the serialization framework.
    """

    cdef CSerialObject* _serial(self):
        return NULL

    cpdef string dumps(
        self,
        str format="asntext",
        bool indent=True,
        bool eol=True,
    ):
        """Dump the object to a byte string.

        Arguments:
            format (`str`): The serialization format to use. Supported formats
                are ``asntext``, ``asnbinary``, ``xml`` and ``json``.
            indent (`bool`): Whether to indent each line in the output.
                Defaults to `True`.
            eol (`bool`): Whether to add newlines in textual formats. 
                ``eol=False`` implies ``indent=False``.

        Returns:
            `bytes`: The serialized object.

        """
        cdef ESerialDataFormat    sdf
        cdef MSerial_Format*      fmt
        cdef string               out
        cdef CNcbiOstream*        s
        cdef CSerialObject*       serial = self._serial()
        cdef stringbuf            buf    = stringbuf()
        cdef TSerial_Format_Flags flags  = 0

        assert serial is not NULL

        if format not in _SERIAL_DATAFORMAT_ENUM:
            raise ValueError(f"invalid format: {format!r}")
        sdf = <ESerialDataFormat> _SERIAL_DATAFORMAT_ENUM[format]

        if not eol:
            indent = False

        if not indent:
            flags |= <TSerial_Format_Flags> ESerial_Xml_Flags.fSerial_Xml_NoIndentation
        if not eol:
            flags |= <TSerial_Format_Flags> ESerial_Xml_Flags.fSerial_Xml_NoEol

        try:
            s = new ostream(&buf)
            fmt = new MSerial_Format(sdf, flags)
            with nogil:
                DumpObject(s[0], serial[0], fmt[0])
            return buf.str()
        finally:
            del s
            del fmt
