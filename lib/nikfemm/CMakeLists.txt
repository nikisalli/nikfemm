
cmake_minimum_required(VERSION 3.0.0)
project(nikfemm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

set (NIKFEMM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set (NIKFEMM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# find all source files in the src directory ending in .cpp or .hpp
file(GLOB_RECURSE NIKFEMM_SRCS ${NIKFEMM_SOURCE_DIR}/*.cpp ${NIKFEMM_SOURCE_DIR}/*.hpp ${NIKFEMM_SOURCE_DIR}/*.h ${NIKFEMM_SOURCE_DIR}/*.c)

include_directories("${NIKFEMM_SOURCE_DIR}")
include_directories("${NIKFEMM_INCLUDE_DIR}")

add_library(nikfemm SHARED ${NIKFEMM_SRCS} ${NIKFEMM_INCLUDE_DIR})

find_package(SDL2 REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${SDL2_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(nikfemm ${SDL2_LIBRARIES})
target_link_libraries(nikfemm ${OpenCV_LIBS})

# get Triangle from lib
# unset TRILIBRARY variable
unset(TRILIBRARY)
add_subdirectory(lib/triangle)
target_link_libraries(nikfemm triangle)

target_include_directories (nikfemm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
