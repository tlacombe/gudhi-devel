cmake_minimum_required(VERSION 2.6)
project(GUDHI)

find_package(Boost		REQUIRED COMPONENTS)

set(CMAKE_CXX_FLAGS "-std=c++0x -g -Wpedantic -Wall")

add_executable(gudh combinatorics/Simplex_tree_node_explicit_storage.cpp main.cpp)
