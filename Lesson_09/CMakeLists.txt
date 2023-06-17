add_executable (Exercise_09.1 Exercise_09.1.cpp Library_09.cpp)
target_link_libraries(Exercise_09.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})