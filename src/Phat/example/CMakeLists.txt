cmake_minimum_required(VERSION 2.6)
project(Phat_interface_examples)

if(PHAT_FOUND)
  add_executable ( rips_persistence_with_phat rips_persistence_with_phat.cpp )
  target_link_libraries(rips_persistence_with_phat ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  add_executable ( cubical_complex_persistence_with_phat cubical_complex_persistence_with_phat.cpp )
  target_link_libraries(cubical_complex_persistence_with_phat ${Boost_SYSTEM_LIBRARY})

  if (TBB_FOUND)
    target_link_libraries(rips_persistence_with_phat ${TBB_LIBRARIES})
    target_link_libraries(cubical_complex_persistence_with_phat ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  add_test(cubical_with_phat_one_sphere ${CMAKE_CURRENT_BINARY_DIR}/cubical_complex_persistence_with_phat ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalOneSphere.txt)
  add_test(cubical_with_phat_two_sphere ${CMAKE_CURRENT_BINARY_DIR}/cubical_complex_persistence_with_phat ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt)

endif(PHAT_FOUND)
