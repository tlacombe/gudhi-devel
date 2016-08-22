cmake_minimum_required(VERSION 2.6)
project(Phat_interface_tests)

if(PHAT_FOUND)

  if (GCOVR_PATH)
    # for gcovr to make coverage reports - Corbera Jenkins plugin
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
  endif()
  if (GPROF_PATH)
    # for gprof to make coverage reports - Jenkins
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
  endif()
  
  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/src/Phat/test/plane_circle" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "${CMAKE_SOURCE_DIR}/src/Phat/test/random_cubical_complex" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_executable ( phat_test phat_test.cpp )
  target_link_libraries(phat_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  
  if (TBB_FOUND)
    target_link_libraries(phat_test ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  # Unitary tests
  add_test(NAME phat_test
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/phat_test 
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/phat_test.xml --log_level=test_suite --report_level=no)

endif(PHAT_FOUND)
