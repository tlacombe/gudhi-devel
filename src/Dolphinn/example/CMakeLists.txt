cmake_minimum_required(VERSION 2.6)
project(Dolphinn_examples)

if (EIGEN3_FOUND)
	add_executable ( Dolphinn_example_radius_from_off Dolphinn_example_radius_from_off.cpp )

	add_executable ( Dolphinn_example_knn_from_points Dolphinn_example_knn_from_points.cpp )

	add_test(NAME Dolphinn_example_radius_from_off COMMAND $<TARGET_FILE:Dolphinn_example_radius_from_off> 
	"${CMAKE_SOURCE_DIR}/data/points/bunny_5000.off" "${CMAKE_SOURCE_DIR}/data/points/bunny_5000.off" "0.1" "0.2" "10")

	add_test(NAME Dolphinn_example_knn_from_points COMMAND $<TARGET_FILE:Dolphinn_example_knn_from_points>
	"20000" "300" "10" "0.2")
endif (EIGEN3_FOUND)
