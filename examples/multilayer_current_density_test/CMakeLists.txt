cmake_minimum_required (VERSION 3.5.0)
project (multilayer_current_density_test CXX)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_executable(multilayer_current_density_test main.cpp)

target_link_libraries (multilayer_current_density_test nikfemm)