enable_testing()

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Add a testing executable
add_executable(unit_tests unit_tests.cpp)

target_link_libraries(unit_tests
    nikfemm
    Boost::unit_test_framework
)

target_compile_definitions(unit_tests
    PRIVATE
        BOOST_TEST_DYN_LINK
)

add_test(test_all unit_tests)