cmake_minimum_required(VERSION 3.20)

if(DEFINED SKBUILD_PROJECT_NAME)
    project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})
else()
    project(pyncbitk-runtime)
endif()

enable_language(C)
enable_language(CXX)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Detect required NCBI C++ Toolkit libraries
#find_package(PCRE REQUIRED)
#find_package(lzo REQUIRED)
#find_package(libunwind REQUIRED)
#find_package(ncbi-cxx-toolkit-public COMPONENTS algo core seqext REQUIRED)

# Configure Toolkit build
set(NCBI_WITHOUT_MT ON)

# Install all dependencies into `pyncbitk-runtime` wheel
if("${SKBUILD_STATE}" STREQUAL "wheel")
    execute_process(
        COMMAND conan install ${CMAKE_SOURCE_DIR} 
            -nr --deployer deployer 
            --deployer-folder "${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME}" 
            -s compiler.cppstd=20)
endif()
