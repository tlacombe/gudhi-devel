cmake_minimum_required(VERSION 2.6)
project(GUDHI)

find_package(Boost		REQUIRED COMPONENTS)

#set(CMAKE_CXX_FLAGS "-std=c++11 -Wpedantic -Wall")

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")


#g++-4.9 -std=c++0x -g -Wpedantic -Wall ../main.cpp;

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

include_directories(. utils combinatorics geometry topology)

add_executable ( TestSkeletonBlockerComplex test/TestSkeletonBlockerComplex.cxx )
add_executable ( TestSimplifiable           test/TestSimplifiable.cxx )
add_executable ( TestContraction            test/TestContraction.cxx )
add_executable ( TestSimplexTree            test/TestSimplexTree.cpp )
add_executable ( TestPersistentCohomology   test/TestPersistentCohomology.cpp )
#add_executable ( demo   test/DemoPersistentCohomology.cpp )
