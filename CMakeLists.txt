cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Numerical_Simulation)

set(CMAKE_CXX_FLAGS "-Wall -O3 --std=c++11")
set(MAKEFLAGS += -s)

#add_subdirectory(general-library)
add_subdirectory(random-library)
add_subdirectory(Lezione_01)
add_subdirectory(Lezione_02)



configure_file( ${CMAKE_SOURCE_DIR}/version_config.h.in ${CMAKE_BINARY_DIR}/version_config.h)