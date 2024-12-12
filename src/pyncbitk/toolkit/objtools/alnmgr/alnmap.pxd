from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.list cimport list as cpplist

from .aln_explorer cimport IAlnExplorer
from ...corelib.ncbiobj cimport CRef, CObject
from ...corelib.ncbimisc cimport TSeqPos, TSignedSeqPos
from ...objects.seqalign.dense_seg cimport CDense_seg, TDim as CDense_seg_TDim, TNumseg as CDense_seg_TNumseg
from ...objects.seqalign.seq_align cimport CSeq_align
from ...objects.seqloc.seq_id cimport CSeq_id


cdef extern from "objtools/alnmgr/alnmap.hpp" namespace "ncbi::objects::CAlnMap" nogil:
    ctypedef unsigned int TSegTypeFlags
    
    enum ESegTypeFlags:
        fSeq                      = 0x0001
        fNotAlignedToSeqOnAnchor  = 0x0002
        fInsert                   = fSeq | fNotAlignedToSeqOnAnchor
        fUnalignedOnRight         = 0x0004
        fUnalignedOnLeft          = 0x0008
        fNoSeqOnRight             = 0x0010
        fNoSeqOnLeft              = 0x0020
        fEndOnRight               = 0x0040
        fEndOnLeft                = 0x0080
        fUnaligned                = 0x0100
        fUnalignedOnRightOnAnchor = 0x0200
        fUnalignedOnLeftOnAnchor  = 0x0400
        fTypeIsSet            
    
    ctypedef CDense_seg_TDim       TDim
    ctypedef TDim                  TNumrow
    # ctypedef CRange[TSeqPos]       TRange
    # ctypedef CRange[TSignedSeqPos] TSignedRange
    ctypedef CDense_seg_TNumseg    TNumseg
    ctypedef cpplist[TSeqPos]      TSeqPosList

    enum EGetChunkFlags:
        fAllChunks           = 0x0000
        fIgnoreUnaligned     = 0x0001
        fInsertSameAsSeq     = 0x0002
        fDeletionSameAsGap   = 0x0004
        fIgnoreAnchor        = fInsertSameAsSeq | fDeletionSameAsGap
        fIgnoreGaps          = 0x0008
        fChunkSameAsSeg      = 0x0010
        fSkipUnalignedGaps   = 0x0020
        fSkipDeletions       = 0x0040
        fSkipAllGaps         = fSkipUnalignedGaps | fSkipDeletions
        fSkipInserts         = 0x0080
        fSkipAlnSeq          = 0x0100
        fSeqOnly             = fSkipAllGaps | fSkipInserts
        fInsertsOnly         = fSkipAllGaps | fSkipAlnSeq
        fAlnSegsOnly         = fSkipInserts | fSkipUnalignedGaps
        fDoNotTruncateSegs   = 0x0200
        fAddUnalignedChunks  = 0x0400

    ctypedef int TGetChunkFlags
    ctypedef TNumseg TNumchunk

    enum ESegmentTrimFlag:
        eSegment_Include
        eSegment_Trim
        eSegment_Remove


cdef extern from "objtools/alnmgr/alnmap.hpp" namespace "ncbi::objects" nogil:

    cppclass CAlnMap(CObject, IAlnExplorer):
        CAlnMap(const CDense_seg& ds)
        CAlnMap(const CDense_seg& ds, TNumrow anchor)

        CRef[CSeq_align] CreateAlignFromRange(
            const vector[TNumrow]& selected_rows,
            TSignedSeqPos          aln_from,
            TSignedSeqPos          aln_to,
            ESegmentTrimFlag       seg_flag = eSegment_Trim)

        # Underlying Dense_seg accessor
        const CDense_seg& GetDenseg() const

        # Dimensions
        TNumseg GetNumSegs() const
        TDim    GetNumRows() const

        # Seq ids
        const CSeq_id& GetSeqId(TNumrow row) const

        # Strands
        bool IsPositiveStrand(TNumrow row) const
        bool IsNegativeStrand(TNumrow row) const
        int  StrandSign      (TNumrow row) const # returns +/- 1

        # Widths
        int  GetWidth        (TNumrow row) const

        # Sequence visible range
        TSignedSeqPos GetSeqAlnStart(TNumrow row) const #aln coords, strand ignored
        TSignedSeqPos GetSeqAlnStop (TNumrow row) const
        # TSignedRange  GetSeqAlnRange(TNumrow row) const
        TSeqPos       GetSeqStart   (TNumrow row) const #seq coords, with strand
        TSeqPos       GetSeqStop    (TNumrow row) const  
        # TRange        GetSeqRange   (TNumrow row) const

        # Segment info
        TSignedSeqPos GetStart  (TNumrow row, TNumseg seg, int offset = 0) const
        TSignedSeqPos GetStop   (TNumrow row, TNumseg seg, int offset = 0) const
        # TSignedRange  GetRange  (TNumrow row, TNumseg seg, int offset = 0) const
        TSeqPos       GetLen    (             TNumseg seg, int offset = 0) const
        TSeqPos       GetSeqLen (TNumrow row, TNumseg seg, int offset = 0) const
        TSegTypeFlags GetSegType(TNumrow row, TNumseg seg, int offset = 0) const
        
        TSegTypeFlags GetTypeAtAlnPos(TNumrow row, TSeqPos aln_pos) const;

        @staticmethod
        bool IsTypeInsert(TSegTypeFlags type)

        # Alignment segments
        TSeqPos GetAlnStart(TNumseg seg) const
        TSeqPos GetAlnStop (TNumseg seg) const
        TSeqPos GetAlnStart()            noexcept const
        TSeqPos GetAlnStop ()            const

        bool    IsSetAnchor()            const
        TNumrow GetAnchor  ()            const
        void    SetAnchor  (TNumrow anchor)
        void    UnsetAnchor()

        # //
        # // Position mapping funcitons
        # // 
        # // Note: Some of the mapping functions have optional parameters
        # //       ESearchDirection dir and bool try_reverse_dir 
        # //       which are used in case an exact match is not found.
        # //       If nothing is found in the ESearchDirection dir and 
        # //       try_reverse_dir == true will search in the opposite dir.

        # TNumseg       GetSeg                 (TSeqPos aln_pos)              const;
        # // if seq_pos falls outside the seq range or into an unaligned region
        # // and dir is provided, will return the first seg in according to dir
        # TNumseg       GetRawSeg              (TNumrow row, TSeqPos seq_pos,
        #                                     ESearchDirection dir = eNone,
        #                                     bool try_reverse_dir = true)  const;
        # // if seq_pos is outside the seq range or within an unaligned region or
        # // within an insert dir/try_reverse_dir will be used
        # TSignedSeqPos GetAlnPosFromSeqPos    (TNumrow row, TSeqPos seq_pos,
        #                                     ESearchDirection dir = eNone,
        #                                     bool try_reverse_dir = true)  const;
        # // if target seq pos is a gap, will use dir/try_reverse_dir
        # TSignedSeqPos GetSeqPosFromSeqPos    (TNumrow for_row,
        #                                     TNumrow row, TSeqPos seq_pos,
        #                                     ESearchDirection dir = eNone,
        #                                     bool try_reverse_dir = true)  const;
        # // if seq pos is a gap, will use dir/try_reverse_dir
        # TSignedSeqPos GetSeqPosFromAlnPos    (TNumrow for_row,
        #                                     TSeqPos aln_pos,
        #                                     ESearchDirection dir = eNone,
        #                                     bool try_reverse_dir = true)  const;
        
        # // Create a vector of relative mapping positions from row0 to row1.
        # // Input:  row0, row1, aln_rng (vertical slice)
        # // Output: result (the resulting vector of positions),
        # //         rng0, rng1 (affected ranges in native sequence coords)
        # void          GetResidueIndexMap     (TNumrow row0,
        #                                     TNumrow row1,
        #                                     TRange aln_rng,
        #                                     vector<TSignedSeqPos>& result,
        #                                     TRange& rng0,
        #                                     TRange& rng1)                 const;