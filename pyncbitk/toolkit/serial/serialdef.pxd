from ..corelib.ncbiobj cimport CObject

cdef extern from "serial/serialdef.hpp" namespace "ncbi" nogil:

    # typedef for object references (constant and nonconstant)
    ctypedef void* TObjectPtr
    ctypedef const void* TConstObjectPtr

    # shortcut typedef: almost everywhere in code we have pointer to const CTypeInfo
    # ctypedef const CTypeInfo* TTypeInfo
    # ctypedef TTypeInfo (*TTypeInfoGetter)(void)
    # ctypedef TTypeInfo (*TTypeInfoGetter1)(TTypeInfo)
    # ctypedef TTypeInfo (*TTypeInfoGetter2)(TTypeInfo, TTypeInfo)
    # ctypedef CTypeInfo* (*TTypeInfoCreator)(void)
    # ctypedef CTypeInfo* (*TTypeInfoCreator1)(TTypeInfo)
    # ctypedef CTypeInfo* (*TTypeInfoCreator2)(TTypeInfo, TTypeInfo)

    # Data file format
    enum ESerialDataFormat: 
        eSerial_None
        eSerial_AsnText  
        eSerial_AsnBinary  
        eSerial_Xml     
        eSerial_Json 

    # Formatting flags
    enum ESerial_AsnText_Flags:
        fSerial_AsnText_NoIndentation
        fSerial_AsnText_NoEol        
    ctypedef unsigned int TSerial_AsnText_Flags

    enum ESerial_Xml_Flags:
        fSerial_Xml_NoIndentation
        fSerial_Xml_NoEol        
        fSerial_Xml_NoXmlDecl     
        fSerial_Xml_NoRefDTD     
        fSerial_Xml_RefSchema   
        fSerial_Xml_NoSchemaLoc  
    ctypedef unsigned int TSerial_Xml_Flags

    enum ESerial_Json_Flags:
        fSerial_Json_NoIndentation
        fSerial_Json_NoEol        
    ctypedef unsigned int TSerial_Json_Flags

    # Data verification parameters
    enum ESerialVerifyData:
        eSerialVerifyData_Default
        eSerialVerifyData_No
        eSerialVerifyData_Never
        eSerialVerifyData_Yes
        eSerialVerifyData_Always
        eSerialVerifyData_DefValue
        eSerialVerifyData_DefValueAlways

    # Skip unknown members parameters
    enum ESerialSkipUnknown:
        eSerialSkipUnknown_Default
        eSerialSkipUnknown_No
        eSerialSkipUnknown_Never
        eSerialSkipUnknown_Yes
        eSerialSkipUnknown_Always

    # File open flags
    enum ESerialOpenFlags:
        eSerial_StdWhenEmpty
        eSerial_StdWhenDash 
        eSerial_StdWhenStd  
        eSerial_StdWhenMask 
        eSerial_StdWhenAny  
        eSerial_UseFileForReread
    ctypedef int TSerialOpenFlags

    # Type family
    enum ETypeFamily:
        eTypeFamilyPrimitive
        eTypeFamilyClass
        eTypeFamilyChoice
        eTypeFamilyContainer
        eTypeFamilyPointer

    # Primitive value type
    enum EPrimitiveValueType:
        ePrimitiveValueSpecial
        ePrimitiveValueBool
        ePrimitiveValueChar
        ePrimitiveValueInteger
        ePrimitiveValueReal
        ePrimitiveValueString
        ePrimitiveValueEnum
        ePrimitiveValueOctetString
        ePrimitiveValueBitString
        ePrimitiveValueAny
        ePrimitiveValueOther

    enum EContainerType:
        eContainerVector
        eContainerList
        eContainerSet
        eContainerMap

    # How to process non-printing character in the ASN VisibleString
    enum EFixNonPrint:
        eFNP_Skip
        eFNP_Allow
        eFNP_Replace
        eFNP_ReplaceAndWarn
        eFNP_Throw
        eFNP_Abort
        eFNP_Default

    # String type
    enum EStringType:
        eStringTypeVisible
        eStringTypeUTF8   

    # How to assign and compare child sub-objects of serial objects
    enum ESerialRecursionMode:
        eRecursive
        eShallow
        eShallowChildless

    # Defines namespace qualification of XML tags
    enum ENsQualifiedMode:
        eNSQNotSet = 0
        eNSUnqualified
        eNSQualified

    enum EDataSpec:
        eUnknown
        eASN    
        eDTD
        eXSD
        eJSON

    enum ESerialFacet:
        eMinLength
        eMaxLength
        eLength
        ePattern

        eInclusiveMinimum
        eExclusiveMinimum
        eInclusiveMaximum
        eExclusiveMaximum
        eMultipleOf

        eMinItems
        eMaxItems
        eUniqueItems

    # Type used for indexing class members and choice variants
    ctypedef size_t TMemberIndex
    ctypedef int TEnumValueType

    # Start if member indexing
    const TMemberIndex kFirstMemberIndex
    # Special value returned from FindMember
    const TMemberIndex kInvalidMember
    # Special value for marking empty choice
    const TMemberIndex kEmptyChoice

    ctypedef ssize_t TPointerOffsetType


    