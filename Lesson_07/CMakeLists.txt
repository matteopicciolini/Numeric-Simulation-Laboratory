add_executable (Exercise_07.1 Exercise_07.1.cpp)
target_link_libraries(Exercise_07.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})