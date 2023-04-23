cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Numerical_Simulation)

set(CMAKE_CXX_FLAGS "-Wall -O3 --std=c++17")
set(MAKEFLAGS += -s)

#add_subdirectory(general-library)
add_subdirectory(random-library)
add_subdirectory(Lesson_01)
add_subdirectory(Lesson_02)
add_subdirectory(Lesson_03)
add_subdirectory(Lesson_04)
add_subdirectory(Lesson_05)
add_subdirectory(Lesson_06)

configure_file( ${CMAKE_SOURCE_DIR}/version_config.h.in ${CMAKE_BINARY_DIR}/version_config.h)