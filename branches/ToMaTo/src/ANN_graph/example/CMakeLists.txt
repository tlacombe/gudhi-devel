cmake_minimum_required(VERSION 2.6)
project(GUDHIAnnGraphExample)

if(ANN_FOUND)
  message("ANN_INCLUDE_DIR = ${ANN_INCLUDE_DIR}")
  include_directories(${ANN_INCLUDE_DIR})

  if (GCOVR_PATH)
    # for gcovr to make coverage reports - Corbera Jenkins plugin
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
  endif()

  add_executable ( simpleanngraph simple_ann_graph.cpp )
  target_link_libraries(simpleanngraph ${ANN_LIBRARIES})
  add_test(simpleanngraph ${CMAKE_CURRENT_BINARY_DIR}/simpleanngraph)
else()
  message("ANN_INCLUDE_DIR NOT FOUND")
endif()
