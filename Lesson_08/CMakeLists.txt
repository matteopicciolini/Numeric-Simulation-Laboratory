add_executable (Exercise_08.1 Exercise_08.1.cpp Library_08.cpp)
target_link_libraries(Exercise_08.1 PRIVATE random-library)

add_executable (Exercise_08.2 Exercise_08.2.cpp Library_08.cpp)
target_link_libraries(Exercise_08.2 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})