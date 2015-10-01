cmake_minimum_required(VERSION 2.6)
project(GUDHIClusterExample)

if(ANN_FOUND)
  message("ANN_INCLUDE_DIR = ${ANN_INCLUDE_DIR}")
  include_directories(${ANN_INCLUDE_DIR})

  if (GCOVR_PATH)
    # for gcovr to make coverage reports - Corbera Jenkins plugin
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
  endif()

  add_executable ( clusteranddensityfromfile cluster_and_density_from_file.cpp )
  target_link_libraries(clusteranddensityfromfile ${ANN_LIBRARIES})
  add_test(clusteranddensityfromfile ${CMAKE_CURRENT_BINARY_DIR}/clusteranddensityfromfile
      ${CMAKE_SOURCE_DIR}/data/points/ToMaTo/toy_example_w_o_density.txt 200 0.25 2 )

  add_executable ( clusterfromfile cluster_from_file.cpp )
  target_link_libraries(clusterfromfile ${ANN_LIBRARIES})
  add_test(clusterfromfile ${CMAKE_CURRENT_BINARY_DIR}/clusterfromfile
      ${CMAKE_SOURCE_DIR}/data/points/ToMaTo/toy_example_w_density.txt 0.25 1e20 )

  add_executable ( anngraphcluster ann_graph_cluster.cpp )
  target_link_libraries(anngraphcluster ${ANN_LIBRARIES})
  add_test(anngraphcluster ${CMAKE_CURRENT_BINARY_DIR}/anngraphcluster)
else()

  message("ANN_INCLUDE_DIR NOT FOUND")
endif()
