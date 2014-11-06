cmake_minimum_required(VERSION 2.6)
project(GUDHI)

find_package(Boost		REQUIRED COMPONENTS)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -Wpedantic -Wall -lboost_system -lgmpxx -lgmp")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

include_directories(. utils combinatorics geometry topology 
                      Nearest_neighbors/include/gudhi
                      Simplex_tree/include/gudhi
                      Persistent_cohomology/include/gudhi 
                      Hasse_complex/include/gudhi 
                      utilities/include/gudhi )
                      
add_subdirectory(test/skbl/)                   


