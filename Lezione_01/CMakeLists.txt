add_executable (Exercise_01.1 Exercise_01.1.cpp Library_01.cpp)
target_link_libraries(Exercise_01.1 PRIVATE random-library)

add_executable (Exercise_01.2 Exercise_01.2.cpp Library_01.cpp)
target_link_libraries(Exercise_01.2 PRIVATE random-library)

add_executable (Exercise_01.3 Exercise_01.3.cpp Library_01.cpp)
target_link_libraries(Exercise_01.3 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories( ${CMAKE_BINARY_DIR})