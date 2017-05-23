cmake_minimum_required(VERSION 2.6)
project(Dolphinn_tests)

include(GUDHI_test_coverage)
if (EIGEN3_FOUND)
	add_executable ( Dolphinn_test Dolphinn_test.cpp )
	target_link_libraries(Dolphinn_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

	gudhi_add_coverage_test(Dolphinn_test)
endif (EIGEN3_FOUND)
