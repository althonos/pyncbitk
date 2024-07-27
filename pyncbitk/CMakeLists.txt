set(CYTHON_DIRECTIVES
    cdivision_warnings
)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

add_custom_command(
  OUTPUT lib.cpp
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/lib.cpp from ${CMAKE_CURRENT_SOURCE_DIR}/example.pyx"
  COMMAND Python::Interpreter -m cython
          "${CMAKE_CURRENT_SOURCE_DIR}/lib.pyx" 
          --output-file lib.cpp
          --cplus
          -X cdivision=True
          -X cdivision_warnings=True
          -X nonecheck=False
          -X warn.undeclared=True
  DEPENDS lib.pyx
  VERBATIM)

python_add_library(lib MODULE WITH_SOABI lib.pyx lib.cpp)

foreach(_lib IN LISTS PYNCBITK_NCBI_LIBS)
    target_link_libraries(lib PUBLIC ${_lib})
endforeach()

install(TARGETS lib DESTINATION ${CMAKE_PROJECT_NAME})

