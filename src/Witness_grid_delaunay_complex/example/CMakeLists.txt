cmake_minimum_required(VERSION 2.6)
SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -frounding-math ")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
include_directories("../include/gudhi")
add_executable ( simpletest simpletest.cpp )
add_executable ( simpletest_alternate simpletest_alternate.cpp )
target_link_libraries( simpletest ${GMP_LIBRARIES} ${CGAL_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries( simpletest_alternate ${GMP_LIBRARIES} ${CGAL_LIBRARY} ${Boost_SYSTEM_LIBRARY})

