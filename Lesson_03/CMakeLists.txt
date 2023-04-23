add_executable (Exercise_03.1 Exercise_03.1.cpp  Library_03.cpp)
target_link_libraries(Exercise_03.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories( ${CMAKE_BINARY_DIR})