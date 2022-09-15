cmake_minimum_required (VERSION 3.0.0)
project (test CXX)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_executable(test main.cpp)

target_link_libraries (test nikfemm)