cmake_minimum_required (VERSION 3.5.0)
project (mtest CXX)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_executable(mtest main.cpp)

target_link_libraries (mtest nikfemm)