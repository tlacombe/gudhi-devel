cmake_minimum_required(VERSION 2.6)
project(Phat_interface_examples)

if(PHAT_FOUND)
  add_executable ( persistence_with_phat_example persistence_with_phat_example.cpp )
  target_link_libraries(persistence_with_phat_example ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  add_test(cubical_with_phat_one_sphere ${CMAKE_CURRENT_BINARY_DIR}/cubical_complex_persistence_with_phat ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalOneSphere.txt)
  add_test(cubical_with_phat_two_sphere ${CMAKE_CURRENT_BINARY_DIR}/cubical_complex_persistence_with_phat ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt)

endif(PHAT_FOUND)
