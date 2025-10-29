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

# Load Python
find_package(Python COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT} REQUIRED)

# Detect the `pyncbitk-runtime` Python library which contains all the build
# dependencies already installed
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import pyncbitk_runtime; print(pyncbitk_runtime.__path__[0], end='')"
    OUTPUT_VARIABLE PYNCBITK_RUNTIME_PATH
)

set(CMAKE_FIND_PACKAGE_REDIRECTS_DIR ${PYNCBITK_RUNTIME_PATH})

# Detect required NCBI C++ Toolkit libraries
find_package(PCRE                    REQUIRED)# PATHS ${PYNCBITK_RUNTIME_PATH} NO_DEFAULT_PATH)
find_package(lzo                     REQUIRED)# PATHS ${PYNCBITK_RUNTIME_PATH} NO_DEFAULT_PATH)
find_package(libunwind               REQUIRED)# PATHS ${PYNCBITK_RUNTIME_PATH} NO_DEFAULT_PATH)
find_package(ncbi-cxx-toolkit-public REQUIRED)# PATHS ${PYNCBITK_RUNTIME_PATH} NO_DEFAULT_PATH)
# find_package(ncbi-cxx-toolkit-public COMPONENTS algo core seqext REQUIRED)

# set(NCBI_PTB_HAS_ROOT ON)
# set(NCBITK_TREE_ROOT ${PYNCBITK_RUNTIME_PATH}/ncbi-cxx-toolkit-public)
# set(NCBITK_SRC_ROOT ${CMAKE_CURRENT_LIST_DIR}/src)
# include(${NCBITK_TREE_ROOT}/res/build-system/cmake/CMake.NCBItoolkit.cmake)

# Configure Toolkit build
# set(NCBI_WITHOUT_MT          ON)
# set(NCBI_PTBCFG_PROJECT_LIST pyncbitk)

# Install shared library objects to `pyncbitk.libs` (like cibuildwheel)
# if(PYNCBITK_INSTALL_LIBS)
#     set(PYTHON_LIBS_INSTALL_DIR "${CMAKE_PROJECT_NAME}.libs")
#     message(STATUS "SKBUILD_PLATLIB_DIR: ${SKBUILD_PLATLIB_DIR}")
#     if("${SKBUILD_STATE}" STREQUAL "wheel")
#         execute_process(COMMAND conan install ${CMAKE_SOURCE_DIR} -nr --deployer runtime_deploy --deployer-folder "${SKBUILD_PLATLIB_DIR}/${PYTHON_LIBS_INSTALL_DIR}")
#     endif()
# endif()

set(PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
include(${CMAKE_CURRENT_LIST_DIR}/src/scripts/cmake/CythonExtension.cmake)

# Set Cython header paths
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include/cython)

# Build project targets
# NCBI_add_subdirectory(src)
add_subdirectory(src)
