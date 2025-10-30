# cython: language_level=3, linetrace=True, binding=True

from libcpp.string cimport string
from libcpp.cast cimport static_cast
from iostream cimport ostream, stringbuf, streambuf

from .toolkit.corelib.ncbistre cimport CNcbiOstream
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


# NOTE: This is a custom `streambuf` that hashes its input using FNV-64 instead 
#       of writing it to a buffer; this allows using the CSerialObject code which 
#       is only implemented for object serialization and use the serialized 
#       representation of the object to hash it on-the-fly.
cdef extern from * nogil:
    """
    class hashbuf : public std::basic_streambuf< char, std::char_traits< char > >
    {
    private:
        class hasher
        {
        private:
            const size_t PRIME = 0x00000100000001b3;  // FNV-64 prime
            const size_t BASIS = 0xcbf29ce484222325;
        public:
            inline hasher() : state(BASIS) {}
            void operator<< (char c) { state ^= c; state *= PRIME; }
            size_t finish() { return state; }
        private:
            size_t state;
        };

        hasher h;

    public:
        inline hashbuf() :
        std::basic_streambuf< char, std::char_traits< char > >(),
        h(hasher())
        {}

        size_t finish() { return h.finish(); }
        std::streamsize xsputn(const char* c, std::streamsize n) override
        {
            for(std::streamsize i = 0; i < n; i++) h << c[i];
            return n;
        }
    };
    """
    cppclass hashbuf(streambuf):
        hashbuf()
        size_t finish() noexcept


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

    def __eq__(self, object other):
        cdef CSerialObject* lostreameft
        cdef CSerialObject* right

        if not isinstance(other, Serial):
            return NotImplemented

        left = self._serial()
        right = (<Serial> other)._serial()

        try:
            return left.Equals(right[0])
        except:
            return NotImplemented

    def __hash__(self):
        cdef ESerialDataFormat    sdf
        cdef MSerial_Format*      fmt
        cdef string               out
        cdef CNcbiOstream*        s
        cdef CSerialObject*       serial = self._serial()
        cdef hashbuf*             hasher
        cdef TSerial_Format_Flags flags  = 0

        assert serial is not NULL

        # use asn binary serialization as it is likely the most compact.
        sdf = ESerialDataFormat.eSerial_AsnBinary
        flags |= <TSerial_Format_Flags> ESerial_Xml_Flags.fSerial_Xml_NoIndentation
        flags |= <TSerial_Format_Flags> ESerial_Xml_Flags.fSerial_Xml_NoEol

        try:
            hasher = new hashbuf()
            s = new ostream(&hasher[0])
            fmt = new MSerial_Format(sdf, flags)
            with nogil:
                DumpObject(s[0], serial[0], fmt[0])
            return hasher.finish()
        finally:
            del s
            del fmt
            del hasher

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
