cmake_minimum_required(VERSION 3.20)

if(DEFINED SKBUILD_PROJECT_NAME)
    project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})
else()
    project(pyncbitk)
endif()

enable_language(C)
enable_language(CXX)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Detect required NCBI C++ Toolkit libraries
find_package(PCRE REQUIRED)
find_package(lzo REQUIRED)
find_package(libunwind REQUIRED)
find_package(ncbi-cxx-toolkit-public COMPONENTS algo core seqext REQUIRED)
NCBIptb_setup()

# Configure Toolkit build
set(NCBI_WITHOUT_MT          ON)
set(NCBI_PTBCFG_PROJECT_LIST pyncbitk)

# Install shared library objects to `pyncbitk.libs` (like cibuildwheel)
if(PYNCBITK_INSTALL_LIBS)
    set(PYTHON_LIBS_INSTALL_DIR "${CMAKE_PROJECT_NAME}.libs")
    message(STATUS "SKBUILD_PLATLIB_DIR: ${SKBUILD_PLATLIB_DIR}")
    if("${SKBUILD_STATE}" STREQUAL "wheel")
        execute_process(COMMAND conan install ${CMAKE_SOURCE_DIR} -nr --deployer runtime_deploy --deployer-folder "${SKBUILD_PLATLIB_DIR}/${PYTHON_LIBS_INSTALL_DIR}")
    endif()
endif()

# Run NCBI C++ Toolkit configuration
include(${CMAKE_CURRENT_LIST_DIR}/src/scripts/cmake/CythonExtension.cmake)

# Set Cython header paths
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include/cython)

# Build project targets
NCBI_add_subdirectory(src)
