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

# Symlink to the C++ Toolkit sources so that NCBIptb is happy with the project structure
file(CREATE_LINK ../vendor/ncbi-cxx-toolkit-public/src     src/ncbi-cxx-toolkit-public     COPY_ON_ERROR SYMBOLIC)
file(CREATE_LINK ../vendor/ncbi-cxx-toolkit-public/include include/ncbi-cxx-toolkit-public COPY_ON_ERROR SYMBOLIC)

# Disable optional components to simplify compilation
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
    LMDB
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
foreach(_component IN LISTS PYNCBITK_DISABLED_COMPONENTS)
    set(NCBI_COMPONENT_${_component}_DISABLED ON)
endforeach()

# Disable SSE4 (only if requested)
if(PYNCBITK_DISABLE_SSE4)
    set(PYNCBITK_PROJECT_FEATURES -OpenMP -SSE)
else()
    set(PYNCBITK_PROJECT_FEATURES -OpenMP)
endif()

# Configure Toolkit build
set(NCBI_WITHOUT_MT          ON)
set(NCBI_PTBCFG_PROJECT_LIST pyncbitk)

# Install shared library objects to `pyncbitk.libs` (like cibuildwheel)
if(PYNCBITK_INSTALL_LIBS)
    set(PYTHON_LIBS_INSTALL_DIR "${CMAKE_PROJECT_NAME}.libs")
endif()

# Run NCBI C++ Toolkit configuration
include(${CMAKE_CURRENT_LIST_DIR}/vendor/ncbi-cxx-toolkit-public/src/build-system/cmake/CMake.NCBItoolkit.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/src/scripts/cmake/CythonExtension.cmake)

# Build project targets
NCBI_add_subdirectory(src)

# Install transitive dependencies to the shared library folders
if(PYNCBITK_INSTALL_LIBS)
    message("Patching install libs")
    get_property(_exts GLOBAL PROPERTY PYNCBITK_CYTHON_EXTENSIONS)
    foreach(_ext IN LISTS _exts)
        get_property(_deps GLOBAL PROPERTY NCBI_PTBPROP_DEPS_${_ext})
        list(APPEND PYNCBITK_DEPENDENCIES ${_deps})
    endforeach()
    list(REMOVE_DUPLICATES PYNCBITK_DEPENDENCIES)
    foreach(_lib IN LISTS PYNCBITK_DEPENDENCIES)
        install(TARGETS ${_lib} DESTINATION ${SKBUILD_PLATLIB_DIR}/${CMAKE_PROJECT_NAME}.libs)
        set_target_properties(${_lib} PROPERTIES INSTALL_RPATH "\$ORIGIN/.")
    endforeach()
endif()
