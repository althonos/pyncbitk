set(CYTHON_DIRECTIVES
    -X cdivision=True
    -X nonecheck=False
)

if(${CMAKE_BUILD_TYPE} EQUAL "Debug")
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
else()
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X boundscheck=True
    -X wraparound=True
  )
endif()

macro(cython_extension _name)
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
            ${CYTHON_DIRECTIVES}
    MAIN_DEPENDENCY 
      ${_name}.pyx
    VERBATIM)
  
  # Add Python library target
  set(_target "py_${_name}")
  python_add_library(${_target} MODULE WITH_SOABI ${_name}.pyx ${_name}.cpp)
  set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${_name} )

  # Link to required C++ Toolkit libraries
  foreach(_lib IN LISTS PYNCBITK_NCBI_LIBS)
    target_link_libraries(${_target} PUBLIC $<TARGET_LINKER_FILE:${_lib}>)
  endforeach()

  # Preserve the relative project structure in the install directory
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS ${_target} DESTINATION ${_dest_folder} )
  message(DEBUG "Install folder for extension ${_name}: ${_dest_folder}")

  # Patch the RPATH to the installed libs (only if libs are installed locally)
  if(DEFINED PYTHON_LIBS_INSTALL_DIR)
    cmake_path(SET _path NORMALIZE ${_dest_folder})
    string(REPLACE "/" ";" _components ${_path})
    set(_rpath "\$ORIGIN/")
    foreach(_x IN LISTS _components)
      string(APPEND _rpath "../")
    endforeach()
    string(APPEND _rpath "${PYTHON_LIBS_INSTALL_DIR}")
    set_target_properties(${_target} PROPERTIES INSTALL_RPATH ${_rpath})
    message(DEBUG "RPATH for extension ${_name}: ${_rpath}")
  endif()

endmacro()