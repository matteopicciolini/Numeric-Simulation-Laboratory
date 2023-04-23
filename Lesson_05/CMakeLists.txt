add_executable (Exercise_05.1 Exercise_05.1.cpp Library_05.cpp)
target_link_libraries(Exercise_05.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})

configure_file(input-output/input.dat input-output/input.dat COPYONLY)