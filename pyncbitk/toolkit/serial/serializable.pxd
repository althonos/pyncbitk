cdef extern from "serial/serializable.hpp" namespace "ncbi::objects::CSerializable" nogil:

    enum EOutputType:
        eAsFasta
        eAsAsnText
        eAsAsnBinary
        eAsXML
        eAsString
    
    cppclass CProxy:
        CProxy(const CSerializable& obj, EOutputType output_type)


cdef extern from "serial/serializable.hpp" namespace "ncbi::objects" nogil:

    cppclass CSerializable:
        CProxy Dump(EOutputType output_type) const



    