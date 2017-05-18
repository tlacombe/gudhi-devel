cmake_minimum_required(VERSION 2.6)
project(Dolphinn_tests)

include(GUDHI_test_coverage)

add_executable ( Dolphinn_test Dolphinn_test.cpp )
target_link_libraries(Dolphinn_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${CGAL_LIBRARY})

gudhi_add_coverage_test(Dolphinn_test)
