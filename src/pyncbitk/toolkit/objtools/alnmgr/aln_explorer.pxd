cdef extern from "objtools/alnmgr/aln_explorer.hpp" namespace "ncbi::IAlnExplorer" nogil:

    enum EAlignType:
        fDNA
        fProtein
        fMixed
        fHomogenous
        fInvalid

    enum ESearchDirection:
        eNone
        eBackwards
        eForward
        eLeft
        eRight

    enum ESortState:
        eUnSorted
        eAscending
        eDescending
        eNotSupported

    # ctypedef CRange[TSeqPos] TRange
    # ctypedef CRange[TSignedSeqPos] TSignedRange


cdef extern from "objtools/alnmgr/aln_explorer.hpp" namespace "ncbi::IAlnSegment" nogil:

    ctypedef unsigned int TSeqTypeFlags

    enum ESeqTypeFlags:
        fAligned
        fGap
        fReversed
        fIndel
        fUnaligned
        fInvalid
        fSeqTypeMask


cdef extern from "objtools/alnmgr/aln_explorer.hpp" namespace "ncbi" nogil:

    cppclass IAlnExplorer:
        pass

    cppclass IAlnSegmentIterator:
        pass
