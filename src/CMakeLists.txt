find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

include("scripts/cmake/CythonExtension.cmake")
set(PYTHON_EXTENSIONS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(pyncbitk)
