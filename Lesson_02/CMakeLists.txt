add_executable (Exercise_02.1 Exercise_02.1.cpp  Library_02.cpp)
target_link_libraries(Exercise_02.1 PRIVATE random-library)

add_executable (Exercise_02.2 Exercise_02.2.cpp  Library_02.cpp)
target_link_libraries(Exercise_02.2 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories( ${CMAKE_BINARY_DIR})