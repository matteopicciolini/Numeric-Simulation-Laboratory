add_executable (Exercise_06.1 Exercise_06.1.cpp)
target_link_libraries(Exercise_06.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})

configure_file(input-output/input.dat input-output/input.dat COPYONLY)
configure_file(run_all_temp.sh run_all_temp.sh COPYONLY)
configure_file(clean.sh clean.sh COPYONLY)