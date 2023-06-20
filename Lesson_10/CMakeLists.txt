find_package(MPI REQUIRED)
add_executable (Exercise_10.1 Exercise_10.1.cpp Library_10.cpp)
target_link_libraries(Exercise_10.1 PRIVATE random-library MPI::MPI_CXX)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})

configure_file(American_capitals.dat American_capitals.dat COPYONLY)