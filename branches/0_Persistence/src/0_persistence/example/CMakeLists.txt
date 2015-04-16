cmake_minimum_required(VERSION 2.6)
project(GUDHI0PersistenceExample)

if(ANN_FOUND)
  message("ANN_INCLUDE_DIR = ${ANN_INCLUDE_DIR}")
  include_directories(${ANN_INCLUDE_DIR})

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
endif()


  add_executable ( simple0persistence simple_0_persistence.cpp )
  target_link_libraries(simple0persistence ${ANN_LIBRARIES})
  add_test(simple0persistence_1 ${CMAKE_CURRENT_BINARY_DIR}/simple0persistence ${CMAKE_SOURCE_DIR}/data/points/0_persistence/spiral_w_density.txt 10 1e-3)
  add_test(simple0persistence_2 ${CMAKE_CURRENT_BINARY_DIR}/simple0persistence ${CMAKE_SOURCE_DIR}/data/points/0_persistence/spiral_w_density.txt 25 1e-3)
  add_test(simple0persistence_0 ${CMAKE_CURRENT_BINARY_DIR}/simple0persistence ${CMAKE_SOURCE_DIR}/data/points/0_persistence/toy_example_w_density.txt 0.25 0)

  add_executable ( densitythen0pers density_then_0_pers.cpp )
  target_link_libraries(densitythen0pers ${ANN_LIBRARIES})
  add_test(densitythen0pers_1 ${CMAKE_CURRENT_BINARY_DIR}/densitythen0pers ${CMAKE_SOURCE_DIR}/data/points/0_persistence/toy_example_w_o_density.txt 200 0.25 1e20)
  add_test(densitythen0pers_2 ${CMAKE_CURRENT_BINARY_DIR}/densitythen0pers ${CMAKE_SOURCE_DIR}/data/points/0_persistence/toy_example_w_o_density.txt 200 0.25 2)
else()
  message("ANN_INCLUDE_DIR NOT FOUND")
endif()
