cmake_minimum_required(VERSION 3.20)

if(DEFINED SKBUILD_PROJECT_NAME)
    project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})
else()
    project(pyncbitk)
endif()

enable_language(C)
enable_language(CXX)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- NCBI C++ Toolkit ---------------------------------------------------------

set(PYNCBITK_PROJECT_LIST
    corelib
    serial
    objects/general
    objects/blastdb
    objects/seq
    objects/seqalign
    objects/seqloc
    objects/seqset
    objtools/readers
    objtools/blast
    objmgr
    algo/blast/api
)

set(PYNCBITK_PROJECT_FEATURES
    -OpenMP
    -SSE
)

set(PYNCBITK_DISABLED_COMPONENTS
    BACKWARD
    BerkeleyDB
    Boost
    Boring
    BZ2
    CASSANDRA
    CURL
    EXSLT
    FASTCGI
    FASTCGIPP
    FONTCONFIG
    FreeType
    FreeTDS
    FTGL
    GCRYPT
    GLEW
    GLPK
    GL2PS
    GMOCK
    GMP
    GPG
    GRPC
    GTK2
    GTEST
    H2O
    KRB5
    JAEGER
    JPEG
    LAPACK
    LEVELDB
    # LMDB
    LZO
    MONGOC
    MONGOCXX
    NCBILS2
    MSGSL
    MySQL
    NCBI
    NCBICRYPT
    NETTLE
    NGHTTP2
    NLohmann
    OpenGL
    OpenSSL
    OPENTRACING
    OSMesa
    PCRE
    Perl
    PerlLibs
    PNG
    PROTOBUF
    PSGLoader
    PubSeqOS
    PYTHON
    ROCKSDB
    SAMTOOLS
    SGE
    # SQLITE3
    Sybase
    TIFF
    THRIFT
    UNWIND
    UUID
    UV
    VDB
    wxWidgets
    WGMLST
    XALAN
    XERCES
    XERCESC
    XML
    XSLT
    XLSXWRITER
    YAML
    Z
    ZSTD
)

set(PYNCBITK_NCBI_LIBS
    blast
    blastdb
    biblio
    # biotree
    bz2
    # cdd
    # cn3d
    connect
    composition_adjustment
    # dbsnp_ptis
    # dbsnp_tooltip_service
    # docsum
    # entrez2
    # entrez2cli
    # entrezgene
    # gbproj
    # gbseq
    # gencoll_client
    gene_info
    gene_info_writer
    general
    # genesbyloc
    genome_collection
    # homologene
    # id1
    # id1cli
    # id2
    # id2cli
    # insdseq
    medline
    # mim
    # m
    # mmdb
    # ncbimime
    # objcoords
    # objprt
    # proj
    pub
    # pubmed
    # regexp
    scoremat
    seq
    seqcode
    seqdb
    seqedit
    seqmasks_io
    seqset
    seqsplit
    ##seqtest
    sequtil
    sqlitewrapp
    submit
    tables
    # tinyseq
    utrtprof
    writedb
    xalgoblastdbindex
    xalgodustmask
    xalgowinmask
    xblast
    xcgi
    xcompress
    xconnect
    xconnext
    # xconnserv
    xcser
    xlogging
    xncbi
    xnetblastcli
    xnetblast
    xobjmgr
    xobjread
    xobjutil
    # xregexp
    xser
    # xthrserv
    xutil
    z
)

# Set path to Toolkit sources
set(NCBITK_TREE_ROOT        ${CMAKE_CURRENT_SOURCE_DIR}/vendor/ncbi-cxx-toolkit-public)
set(NCBITK_SRC_ROOT         ${CMAKE_CURRENT_SOURCE_DIR}/vendor/ncbi-cxx-toolkit-public/src)
set(NCBI_TREE_ROOT          ${CMAKE_CURRENT_SOURCE_DIR}/vendor/ncbi-cxx-toolkit-public)
set(NCBI_SRC_ROOT           ${CMAKE_CURRENT_SOURCE_DIR}/vendor/ncbi-cxx-toolkit-public/src)
set(NCBI_BUILD_ROOT         ${CMAKE_CURRENT_BINARY_DIR}/ncbi)

# Configure Toolkit build
set(NCBI_WITHOUT_MT                ON)
set(NCBI_VERBOSE_PROJECTS          "")
set(NCBI_PTBCFG_PROJECT_COMPONENTS "")
set(NCBI_PTBCFG_PROJECT_FEATURES   ${PYNCBITK_PROJECT_FEATURES})
set(NCBI_PTBCFG_PROJECT_LIST       ${PYNCBITK_PROJECT_LIST})
set(NCBI_PTBCFG_PROJECT_TAGS       "-test")
set(NCBI_PTBCFG_PROJECT_TARGETS    "")
set(NCBI_PTBCFG_ALLOW_COMPOSITE    OFF)
set(NCBI_PTBCFG_SKIP_ANALYSIS      OFF)

# Disable unneeded external dependencies
foreach(_component IN LISTS PYNCBITK_DISABLED_COMPONENTS)
    set(NCBI_COMPONENT_${_component}_DISABLED ON)
endforeach()

# Run NCBI C++ Toolkit configuration and build libraries
include(${CMAKE_CURRENT_SOURCE_DIR}/vendor/ncbi-cxx-toolkit-public/src/build-system/cmake/CMake.NCBItoolkit.cmake)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/ncbi-cxx-toolkit-public/src)

# Register library directory
link_directories(${CMAKE_BINARY_DIR}/../lib)

# Install required libraries
if(PYNCBITK_INSTALL_LIBS)
    set(PYTHON_LIBS_INSTALL_DIR "${CMAKE_PROJECT_NAME}.libs")
    foreach(_lib IN LISTS PYNCBITK_NCBI_LIBS)
        if(TARGET ${_lib})
            install(TARGETS ${_lib} DESTINATION ${SKBUILD_PLATLIB_DIR}/${CMAKE_PROJECT_NAME}.libs)
        else()
            message(SEND_ERROR "Unknown requested library: ${_lib}")
        endif()
        set_target_properties(${_lib} PROPERTIES INSTALL_RPATH "\$ORIGIN/.")
    endforeach()
endif()

# --- PyNCBItk -----------------------------------------------------------------

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Build Python extension
add_subdirectory(src)
