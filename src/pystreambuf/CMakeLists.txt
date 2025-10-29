add_library(pystreambuf STATIC pyreadbuf.cpp pyreadintobuf.cpp pywritebuf.cpp)
target_include_directories(pystreambuf PUBLIC ${Python_INCLUDE_DIRS} ${CMAKE_CURRENT_LIST_DIR})
target_link_directories(pystreambuf PUBLIC ${CMAKE_CURRENT_BINARY_DIR})