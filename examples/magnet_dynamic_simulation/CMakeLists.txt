cmake_minimum_required (VERSION 3.0.0)
project (magnet_dynamic_simulation CXX)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_executable(magnet_dynamic_simulation main.cpp)

target_link_libraries (magnet_dynamic_simulation nikfemm)