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
set(CMAKE_CXX_STANDARD 20)

# Load Python
find_package(Python COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT} REQUIRED)

# Detect the `pyncbitk-runtime` Python library which contains all the build
# dependencies already installed
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import pyncbitk_runtime; print(pyncbitk_runtime.__path__[0], end='')"
    OUTPUT_VARIABLE PYNCBITK_RUNTIME_PATH
    COMMAND_ERROR_IS_FATAL ANY
)
cmake_path(GET PYNCBITK_RUNTIME_PATH PARENT_PATH SITE_PACKAGES_PATH)

# Detect required libraries only in the given folder
set(CMAKE_FIND_PACKAGE_REDIRECTS_DIR ${PYNCBITK_RUNTIME_PATH})
if(NOT APPLE)
    find_package(libunwind REQUIRED)
endif()
find_package(PCRE REQUIRED)
find_package(lzo REQUIRED)
find_package(ncbi-cxx-toolkit-public REQUIRED)

# Collect the runtime library folders (to use to set $RPATH)
set(RUNTIME_LIBRARY_DIRS)
foreach(_lib IN LISTS ncbi-cxx-toolkit-public_LIBRARIES)
  string(REGEX REPLACE "^${SITE_PACKAGES_PATH}/?" "" _lib_path ${_lib})
  cmake_path(GET _lib_path PARENT_PATH _lib_folder)
  list(APPEND RUNTIME_LIBRARY_DIRS ${_lib_folder})
endforeach()
list(REMOVE_DUPLICATES RUNTIME_LIBRARY_DIRS)

# Set path to the Python extension directory
set(PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
include(${CMAKE_CURRENT_LIST_DIR}/src/scripts/cmake/CythonExtension.cmake)

# Set path to the  Cython headers
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include/cython)

# Build project targets
add_subdirectory(src)
