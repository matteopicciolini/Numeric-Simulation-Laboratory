add_executable (Exercise_04.1 Exercise_04.1.cpp)
target_link_libraries(Exercise_04.1 PRIVATE random-library)

add_definitions(-D_VAR=\"${CMAKE_SOURCE_DIR}\")
include_directories(${CMAKE_BINARY_DIR})

configure_file(input/eq/gas/input.gas input/eq/gas/input.gas COPYONLY)
configure_file(input/eq/liquid/input.liquid input/eq/liquid/input.liquid COPYONLY)
configure_file(input/eq/solid/input.solid input/eq/solid/input.solid COPYONLY)

configure_file(input/input.gas input/input.gas COPYONLY)
configure_file(input/input.liquid input/input.liquid COPYONLY)
configure_file(input/input.solid input/input.solid COPYONLY)

configure_file(input/eq/config.in input/eq/config.in COPYONLY)