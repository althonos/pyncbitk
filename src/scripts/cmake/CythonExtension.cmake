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
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${NCBI_CURRENT_SOURCE_DIR})
  string(REPLACE "/" "." _target ${_dest_folder}.${_name})

  # Build the Python extension as an NCBIptb custom target
  function(${_name}_definition)
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
    target_link_libraries(${_target} PUBLIC ${NCBITMP_NCBILIB} ${NCBITMP_EXTLIB})
    foreach(_dep IN LISTS CYTHON_EXTENSION_DEPENDS)
      if(TARGET ${_dep})
        target_include_directories(${_target} PUBLIC $<TARGET_PROPERTY:${_dep},INCLUDE_DIRECTORIES>)
      endif()
    endforeach()

    # Preserve the relative project structure in the install directory
    string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
    install(TARGETS ${_target} DESTINATION ${_dest_folder} )
    message(DEBUG "Install folder for extension ${_name}: ${_dest_folder}")

    # Patch the RPATH to the installed libs (only if libs are installed locally)
    if(PYNCBITK_INSTALL_LIBS AND DEFINED PYTHON_LIBS_INSTALL_DIR)
      cmake_path(SET _path NORMALIZE ${_dest_folder})
      string(REPLACE "/" ";" _components ${_path})
      set(_rpath "\$ORIGIN/")
      foreach(_x IN LISTS _components)
        string(APPEND _rpath "../")
      endforeach()
      string(APPEND _rpath "${PYTHON_LIBS_INSTALL_DIR}")
      set_target_properties(${_target} PROPERTIES INSTALL_RPATH ${_rpath})
      message(DEBUG "RPATH for extension ${_name}: ${_rpath}")
    else()
      set_target_properties(${_target} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif()

  endfunction()
  NCBI_begin_custom_target(${_target})
    NCBI_project_tags(python)
    NCBI_custom_target_definition(${_name}_definition)
    NCBI_uses_toolkit_libraries(${CYTHON_EXTENSION_DEPENDS})
  NCBI_end_custom_target()

  # Add the targets to the list of Cython extensions
  get_property(_ext GLOBAL PROPERTY PYNCBITK_CYTHON_EXTENSIONS)
  list(APPEND _ext ${_target})
  set_property(GLOBAL PROPERTY PYNCBITK_CYTHON_EXTENSIONS ${_ext})
endmacro()
