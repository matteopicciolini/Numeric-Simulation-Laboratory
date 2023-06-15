add_executable (Exercise_07.1 Exercise_07.1.cpp)
target_link_libraries(Exercise_07.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})

configure_file(input-output/eq/gas/input.gas input-output/eq/gas/input.gas COPYONLY)
configure_file(input-output/eq/liquid/input.liquid input-output/eq/liquid/input.liquid COPYONLY)
configure_file(input-output/eq/solid/input.solid input-output/eq/solid/input.solid COPYONLY)

configure_file(input-output/input.gas input-output/input.gas COPYONLY)
configure_file(input-output/input.liquid input-output/input.liquid COPYONLY)
configure_file(input-output/input.solid input-output/input.solid COPYONLY)

configure_file(input-output/eq/config.in input-output/eq/config.in COPYONLY)