cmake_minimum_required (VERSION 3.0.0)
project (gradient_descent_design_optimization CXX)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_executable(gradient_descent_design_optimization main.cpp)

target_link_libraries (gradient_descent_design_optimization nikfemm)