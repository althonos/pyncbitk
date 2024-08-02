from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..corelib.ncbiobj cimport CObject, CRef

cdef extern from "corelib/ncbi_tree.hpp" namespace "ncbi" nogil:

    cppclass CTreeNode[TValue, TKeyGetterP]:
        pass

    cppclass CPairNodeKeyGetter[TNode, TKeyEquapP]:
        ctypedef TNode TNodeTyp

    cppclass CTreePair[TId, TValue, TIdEqual]:
        ctypedef TId     TIdType
        ctypedef TValue  TValueType

        ctypedef CTreePair[TId, TValue, TIdEqual]        TTreePair
        ctypedef CPairNodeKeyGetter[TTreePair, TIdEqual] TPairKeyGetter
        ctypedef CTreeNode[TTreePair, TPairKeyGetter]    TPairTreeNode