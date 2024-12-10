from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "objtools/align_format/format_flags.hpp" namespace "ncbi::align_format" nogil:

    const string kArgOutputFormat
    const int kDfltArgOutputFormat
    const string kArgShowGIs
    const bool kDfltArgShowGIs
    const string kArgNumDescriptions
    const size_t kDfltArgNumDescriptions
    const string kArgNumAlignments
    const size_t kDfltArgNumAlignments
    const string kArgProduceHtml
    const bool kDfltArgProduceHtml