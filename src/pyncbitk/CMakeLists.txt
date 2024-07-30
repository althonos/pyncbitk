find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# --- Cython compilation -------------------------------------------------------

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
  python_add_library(${_name} MODULE WITH_SOABI ${_name}.pyx ${_name}.cpp)
  foreach(_lib IN LISTS PYNCBITK_NCBI_LIBS)
    target_link_libraries(${_name} PUBLIC ${_lib})
  endforeach()
  install(TARGETS ${_name} DESTINATION ${CMAKE_PROJECT_NAME})
endmacro()

cython_extension(objects)
cython_extension(objmgr)
cython_extension(objtools)
cython_extension(algo)
