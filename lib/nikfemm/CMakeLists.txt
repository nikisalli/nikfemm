
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
add_library(nikfemm_static STATIC ${NIKFEMM_SRCS} ${NIKFEMM_INCLUDE_DIR})

if(NOT DEFINED NIKFEMM_USE_OPENCV)
    set(NIKFEMM_USE_OPENCV ON)
endif()

# add opencv if NIKFEMM_USE_OPENCV option is set
if (NIKFEMM_USE_OPENCV)
    find_package(OpenCV REQUIRED)
    target_link_libraries(nikfemm ${OpenCV_LIBS})
    target_compile_definitions(nikfemm PUBLIC NIKFEMM_USE_OPENCV)
endif()

target_include_directories (nikfemm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# add tests for make test
if (NOT DEFINED NIKFEMM_BUILD_TESTS)
    set(NIKFEMM_BUILD_TESTS ON)
endif()

if (NIKFEMM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
