
cmake_minimum_required(VERSION 3.5.0)
project(nikfemm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -mavx")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (NIKFEMM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set (NIKFEMM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# find all source files in the src directory ending in .cpp or .hpp
file(GLOB_RECURSE NIKFEMM_SRCS ${NIKFEMM_SOURCE_DIR}/*.cpp ${NIKFEMM_SOURCE_DIR}/*.hpp ${NIKFEMM_SOURCE_DIR}/*.h ${NIKFEMM_SOURCE_DIR}/*.c)

include_directories("${NIKFEMM_SOURCE_DIR}")
include_directories("${NIKFEMM_INCLUDE_DIR}")

add_library(nikfemm SHARED ${NIKFEMM_SRCS} ${NIKFEMM_INCLUDE_DIR})

# add opencv if OPENCV option is set
option(OPENCV "Use OpenCV" OFF)
if (OPENCV)
    add_definitions(-DNIKFEMM_USE_OPENCV)
    find_package(OpenCV REQUIRED)
    include_directories(${OpenCV_INCLUDE_DIRS})
    target_link_libraries(nikfemm ${OpenCV_LIBS})
endif()

target_include_directories (nikfemm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# tests
add_subdirectory(test)
