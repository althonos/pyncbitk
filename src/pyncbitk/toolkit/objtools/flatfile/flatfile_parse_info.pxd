from libcpp cimport bool


cdef extern from "objtools/flatfile/flatfile_parse_info.hpp" namespace "ncbi::Parser::EOutput" nogil:

    enum EOutput:
        BioSeqSet
        Seqsubmit


cdef extern from "objtools/flatfile/flatfile_parse_info.hpp" namespace "ncbi::Parser::EMode" nogil:

    enum EMode:
        Release
        HTGS
        HTGSCON
        Relaxed


cdef extern from "objtools/flatfile/flatfile_parse_info.hpp" namespace "ncbi::Parser::ESource" nogil:

    enum ESource:
        unknown
        NCBI
        EMBL
        GenBank
        DDBJ
        LANL
        SPROT
        RefSeq
        Flybase
        USPTO
        All


cdef extern from "objtools/flatfile/flatfile_parse_info.hpp" namespace "ncbi::Parser::EFormat" nogil:

    enum EFormat:
        unknown
        EMBL
        GenBank
        SPROT
        DDBJ
        XML
        ALL


cdef extern from "objtools/flatfile/flatfile_parse_info.hpp" namespace "ncbi" nogil:

    cppclass FileBuf:
        FileBuf()

    cppclass Parser:
        Parser()

        EOutput output_format

        EFormat format
        ESource source
        bool    all

        EMode   mode

        void InitializeKeywordParser(EFormat) except +
        # CKeywordParser& KeywordParser();