
cdef extern from "algo/blast/blastinput/blast_args.hpp" namespace "ncbi::blast::CFormattingArgs" nogil:

    enum EOutputFormat:
        ePairwise
        eQueryAnchoredIdentities
        eQueryAnchoredNoIdentities
        eFlatQueryAnchoredIdentities
        eFlatQueryAnchoredNoIdentities
        eXml
        eTabular
        eTabularWithComments
        eAsnText
        eAsnBinary
        eCommaSeparatedValues
        eArchiveFormat
        eJsonSeqalign
        eJson
        eXml2
        eJson_S
        eXml2_S
        eSAM
        eTaxFormat
        eAirrRearrangement
        eFasta
        eEndValue

    enum EFormatFlags:
        eDefaultFlag
        eIsVDB
        eIsSAM
        eIsVDB_SAM
        eIsAirrRearrangement