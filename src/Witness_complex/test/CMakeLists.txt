cmake_minimum_required(VERSION 2.6)
project(Witness_complex_tests)

include(GUDHI_test_coverage)

add_executable ( Witness_complex_test_simple_witness_complex test_simple_witness_complex.cpp )
target_link_libraries(Witness_complex_test_simple_witness_complex ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Witness_complex_test_simple_witness_complex ${TBB_LIBRARIES})
endif(TBB_FOUND)

gudhi_add_coverage_test(Witness_complex_test_simple_witness_complex)

add_executable ( Witness_complex_test_implementation_comparison test_implementation_comparison.cpp )
target_link_libraries(Witness_complex_test_implementation_comparison ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(Witness_complex_test_implementation_comparison ${TBB_LIBRARIES})
# CGAL and Eigen3 are required for Euclidean version of Witness
if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( Witness_complex_test_euclidean_simple_witness_complex test_euclidean_simple_witness_complex.cpp )
  target_link_libraries(Witness_complex_test_euclidean_simple_witness_complex ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Witness_complex_test_euclidean_simple_witness_complex ${TBB_LIBRARIES})
  endif(TBB_FOUND)

add_executable ( generate_landmarks generate_landmarks.cpp )
target_link_libraries(generate_landmarks ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(generate_landmarks ${TBB_LIBRARIES})

add_executable ( generate_sphere generate_sphere.cpp )
target_link_libraries(generate_sphere ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(generate_sphere ${TBB_LIBRARIES})
      # Unitary tests definition and xml result file generation
      add_test(NAME euclidean_simple_witness_complex
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Witness_complex_test_euclidean_simple_witness_complex
          # XML format for Jenkins xUnit plugin
          --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Witness_complex_test_euclidean_simple_witness_complexUT.xml --log_level=test_suite --report_level=no)
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.6.0)
endif()
  gudhi_add_coverage_test(Witness_complex_test_euclidean_simple_witness_complex)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
