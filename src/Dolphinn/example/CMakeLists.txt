cmake_minimum_required(VERSION 2.6)
project(Dolphinn_examples)

	add_executable ( Dolphinn_example_radius_from_off Dolphinn_example_radius_from_off.cpp )
	target_link_libraries(Dolphinn_example_radius_from_off  ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${CGAL_LIBRARY})

	add_executable ( Dolphinn_example_knn_from_points Dolphinn_example_knn_from_points.cpp )
	target_link_libraries(Dolphinn_example_knn_from_points  ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${CGAL_LIBRARY})

	add_test(NAME Dolphinn_example_radius_from_off COMMAND $<TARGET_FILE:Dolphinn_example_radius_from_off> 
	"${CMAKE_SOURCE_DIR}/data/points/bunny_5000.off" "${CMAKE_SOURCE_DIR}/data/points/bunny_5000.off" "0.1" "0.2" "10")

	add_test(NAME Dolphinn_example_knn_from_points COMMAND $<TARGET_FILE:Dolphinn_example_knn_from_points>
	"20000" "300" "10")
