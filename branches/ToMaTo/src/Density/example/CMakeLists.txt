cmake_minimum_required(VERSION 2.6)
project(GUDHIDensityExample)

if(ANN_FOUND)
  message("ANN_INCLUDE_DIR = ${ANN_INCLUDE_DIR}")
  include_directories(${ANN_INCLUDE_DIR})

  if (GCOVR_PATH)
    # for gcovr to make coverage reports - Corbera Jenkins plugin
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
  endif()

  add_executable ( anngraphdensity ann_graph_density.cpp )
  target_link_libraries(anngraphdensity ${ANN_LIBRARIES})
  add_test(anngraphdensity ${CMAKE_CURRENT_BINARY_DIR}/anngraphdensity)
else()
  message("ANN_INCLUDE_DIR NOT FOUND")
endif()
