cmake_minimum_required(VERSION 2.6)
project(Topological_inference_with_cubical_complexes_test)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( Topological_inference_with_cubical_complexesCCUT Topological_inference_with_cubical_complexes_test.cpp )
target_link_libraries(Topological_inference_with_cubical_complexesCCUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Topological_inference_with_cubical_complexesCCUT ${TBB_LIBRARIES})
endif()

#just to test Clement's kd trees code: periodic_kd_tree_test.cpp. Remove when it is confirmed that this work well
#add_executable ( periodic_kd_tree_testCCUT periodic_kd_tree_test.cpp )
#target_link_libraries(periodic_kd_tree_testCCUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#if (TBB_FOUND)
#  target_link_libraries(periodic_kd_tree_testCCUT ${TBB_LIBRARIES})
#endif()


# Unitary tests
add_test(NAME Topological_inference_with_cubical_complexesCCUT
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Topological_inference_with_cubical_complexesCCUT 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Topological_inference_with_cubical_complexesCCUT.xml --log_level=test_suite --report_level=no)

