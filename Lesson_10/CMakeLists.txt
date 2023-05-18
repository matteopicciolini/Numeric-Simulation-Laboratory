add_executable (Exercise_10.1 Exercise_10.1.cpp)
target_link_libraries(Exercise_10.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})