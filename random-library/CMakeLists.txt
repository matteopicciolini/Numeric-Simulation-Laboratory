add_library(
    random-library
    random.h
    random.cpp
)
target_include_directories(random-library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})