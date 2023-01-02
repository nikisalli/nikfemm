cmake_minimum_required (VERSION 3.0.0)
project (integration_test CXX)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_executable(integration_test main.cpp)

target_link_libraries (integration_test nikfemm)