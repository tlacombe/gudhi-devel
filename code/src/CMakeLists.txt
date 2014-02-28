cmake_minimum_required(VERSION 2.6)
project(GUDHI)

find_package(Boost		REQUIRED COMPONENTS)

set(CMAKE_CXX_FLAGS "-std=c++0x -g -Wpedantic -Wall")

#g++-4.9 -std=c++0x -g -Wpedantic -Wall ../main.cpp;


include_directories(. utils combinatorics geometry)



#add_executable(
#	TestSkeletonBlockerComplex test/TestSkeletonBlockerComplex.cxx
#)

#add_executable(
#	TestSimplifiable test/TestSimplifiable.cxx
#)

#add_executable(
#	TestContraction test/TestContraction.cxx
#)

add_executable(TestSimplexTree test_simplex_tree.cpp)
