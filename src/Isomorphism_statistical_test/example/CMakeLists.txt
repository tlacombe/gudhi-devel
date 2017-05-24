cmake_minimum_required(VERSION 2.6)
project(Isomorphism_statistical_examples)

add_executable (example_read_distance_matrix example_read_distance_matrix.cpp)
add_executable (example_simple example_simple.cpp)
add_executable (example_step_by_step example_step_by_step.cpp)

target_link_libraries(example_read_distance_matrix ${Boost_LIBRARIES})
target_link_libraries(example_simple ${Boost_LIBRARIES})
target_link_libraries(example_step_by_step ${Boost_LIBRARIES})
