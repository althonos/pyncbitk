find_package(Python COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT} REQUIRED)
set(PYTHON_EXTENSIONS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

# --- Detect PyInterpreterState_GetID ------------------------------------------

include(CheckSymbolExists)

set(SAFE_CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(SAFE_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
set(SAFE_CMAKE_REQUIRED_LINK_DIRECTORIES "${CMAKE_REQUIRED_LINK_DIRECTORIES}")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${Python_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${Python_LIBRARIES})
set(CMAKE_REQUIRED_LINK_DIRECTORIES ${CMAKE_REQUIRED_LINK_DIRECTORIES} ${Python_LIBRARY_DIRS})
check_symbol_exists(PyInterpreterState_GetID "stdint.h;stdlib.h;Python.h" HAVE_PYINTERPRETERSTATE_GETID)
set(CMAKE_REQUIRED_INCLUDES "${SAFE_CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_LIBRARIES "${SAFE_CMAKE_REQUIRED_LIBRARIES}")
set(CMAKE_REQUIRED_LINK_DIRECTORIES "${SAFE_CMAKE_REQUIRED_LINK_DIRECTORIES}")
set(PYSTATE_PATCH_H ${CMAKE_CURRENT_LIST_DIR}/pystate_patch.h)

# --- Prepare Cython directives and constants ----------------------------------

set(CYTHON_DIRECTIVES
    -X cdivision=True
    -X nonecheck=False
    -E SYS_IMPLEMENTATION_NAME=$<LOWER_CASE:${Python_INTERPRETER_ID}>
    -E SYS_VERSION_INFO_MAJOR=${Python_VERSION_MAJOR}
    -E SYS_VERSION_INFO_MINOR=${Python_VERSION_MINOR}
    -E TARGET_CPU=$<LOWER_CASE:${CMAKE_SYSTEM_PROCESSOR}>
    -E PROJECT_VERSION=${CMAKE_PROJECT_VERSION}
)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X cdivision_warnings=True
    -X warn.undeclared=True
    -X warn.unreachable=True
    -X warn.maybe_uninitialized=True
    -X warn.unused=True
    -X warn.unused_arg=True
    -X warn.unused_result=True
    -X warn.multiple_declarators=True
  )
  if(NOT Python_INTERPRETER_ID STREQUAL PyPy)
    set(CYTHON_DIRECTIVES
      ${CYTHON_DIRECTIVES}
      -X linetrace=true
    )
  endif()
else()
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X boundscheck=False
    -X wraparound=False
  )
endif()

if((NOT "${SKBUILD_SABI_VERSION}" STREQUAL "") AND (NOT CMAKE_BUILD_TYPE STREQUAL Debug) AND (NOT SKBUILD_STATE STREQUAL editable))
  message(STATUS "Building in Limited API mode for Python: ${SKBUILD_SABI_VERSION}")
else()
  message(STATUS "Building in latest API mode for Python: ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}")
endif()

macro(cython_extension _name)
  set(multiValueArgs DEPENDS LINKS)
  cmake_parse_arguments(CYTHON_EXTENSION "" "" "${multiValueArgs}" ${ARGN} )

  # Make sure that the source directory is known
  if(NOT DEFINED PYTHON_EXTENSIONS_SOURCE_DIR)
    message(FATAL_ERROR "The PYTHON_EXTENSIONS_SOURCE_DIR variable has not been set.")
  endif()

  # Generate C++ file from Cython file
  add_custom_command(
    OUTPUT ${_name}.cpp
    COMMENT
      "Making ${CMAKE_CURRENT_BINARY_DIR}/${_name}.cpp from ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.pyx"
    COMMAND
      Python::Interpreter -m cython
            "${CMAKE_CURRENT_SOURCE_DIR}/${_name}.pyx"
            --output-file ${_name}.cpp
            --cplus
            --depfile
            -I "${CYTHON_HEADERS_DIR}"
            ${CYTHON_DIRECTIVES}
    MAIN_DEPENDENCY
      ${_name}.pyx
    DEPFILE
      ${_name}.cpp.dep
    VERBATIM)

  # Build fully-qualified module name as the target name
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_LIST_DIR})
  string(REPLACE "/" "." _target ${_dest_folder}.${_name})

  # Add Python library target
  if((NOT "${SKBUILD_SABI_VERSION}" STREQUAL "") AND (NOT CMAKE_BUILD_TYPE STREQUAL Debug) AND (NOT SKBUILD_STATE STREQUAL editable))
    python_add_library(${_target} MODULE WITH_SOABI USE_SABI "${SKBUILD_SABI_VERSION}" ${_name}.pyx ${_name}.cpp)
  else()
    python_add_library(${_target} MODULE WITH_SOABI ${_name}.pyx ${_name}.cpp)
  endif()
  set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${_name} )

  # Add debug flags
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    if(NOT Python_INTERPRETER_ID STREQUAL PyPy)
      target_compile_definitions(${_target} PUBLIC CYTHON_TRACE_NOGIL=1)
    endif()
  else()
    target_compile_definitions(${_target} PUBLIC CYTHON_WITHOUT_ASSERTIONS=1)
  endif()

  # Include patch for `PyInterpreterState_GetID` to all Python extensions
  target_precompile_headers(${_target} PRIVATE ${PYSTATE_PATCH_H})

  # Link to NCBI libraries and add include directories if needed
  # message(STATUS "libs NCBI=(${NCBITMP_NCBILIB}) EXT=(${NCBITMP_EXTLIB})")
  # message(STATUS "${ncbi-cxx-toolkit-public_INCLUDE_DIRS}")
  # message(STATUS "${ncbi-cxx-toolkit-public_LIBRARIES}")

  # target_link_libraries(${_target} PUBLIC ${NCBITMP_NCBILIB} ${NCBITMP_EXTLIB})
  target_link_libraries(${_target} PUBLIC pystreambuf)

  target_include_directories(${_target} PUBLIC ${ncbi-cxx-toolkit-public_INCLUDE_DIRS})
  target_link_libraries(${_target} PUBLIC ${ncbi-cxx-toolkit-public_LIBRARIES})

  # target_link_directories(${_target} PUBLIC "/home/althonos/.local/lib/python3.13/site-packages/pyncbitk_runtime/ncbi-cxx-toolkit-public/lib" )
  foreach(_dep IN LISTS CYTHON_EXTENSION_DEPENDS)
    # message(STATUS "dep: ${_dep}")
    # NCBI_internal_identify_libs(_link _dep2)
    # message(STATUS "link: ${_link} dep2: ${_dep2} dep: ${_dep}")
    target_link_libraries(${_target} PUBLIC ${_dep})
    if(TARGET ${_dep})
      target_include_directories(${_target} PUBLIC $<TARGET_PROPERTY:${_dep},INCLUDE_DIRECTORIES>)
    endif()
  endforeach()

  # Preserve the relative project structure in the install directory
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS ${_target} DESTINATION ${_dest_folder} )
  message(DEBUG "Install folder for extension ${_name}: ${_dest_folder}")
  message(DEBUG "(${_target}) setting install folder: ${_dest_folder}")

  # Patch the RPATH to the installed libs using relative paths 
  cmake_path(SET _path NORMALIZE ${_dest_folder})
  string(REPLACE "/" ";" _components ${_path})
  if(APPLE)
    set(_rpath "@loader_path/")
  else()
    set(_rpath "\$ORIGIN/")
  endif()
  foreach(_x IN LISTS _components)
    string(APPEND _rpath "../")
  endforeach()

  # ensure all folders are added (if multiple folders)
  foreach(_folder IN LISTS RUNTIME_LIBRARY_DIRS)
    set_target_properties(${_target} PROPERTIES INSTALL_RPATH "${_rpath}${_folder}")
  endforeach()

  # Add the targets to the list of Cython extensions
  get_property(_ext GLOBAL PROPERTY PYNCBITK_CYTHON_EXTENSIONS)
  list(APPEND _ext ${_target})
  set_property(GLOBAL PROPERTY PYNCBITK_CYTHON_EXTENSIONS ${_ext})
endmacro()
