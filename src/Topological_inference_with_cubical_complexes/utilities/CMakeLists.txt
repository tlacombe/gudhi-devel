cmake_minimum_required(VERSION 2.6)
project(GUDHITopologicalInferenceWithCubicalComplexes)

add_executable ( topological_inference topological_inference.cpp )
target_link_libraries(topological_inference ${Boost_SYSTEM_LIBRARY})

add_executable ( topological_inference_periodic_domain topological_inference_periodic_domain.cpp )
target_link_libraries(topological_inference_periodic_domain ${Boost_SYSTEM_LIBRARY})
   
add_executable ( morphological_operation_on_cubical_complex morphological_operation_on_cubical_complex.cpp )	
target_link_libraries(morphological_operation_on_cubical_complex ${Boost_SYSTEM_LIBRARY})

add_executable ( morphological_operation_distance_from_point_cloud morphological_operation_distance_from_point_cloud.cpp )
target_link_libraries(morphological_operation_distance_from_point_cloud ${Boost_SYSTEM_LIBRARY})


if (CGAL_FOUND)
	target_link_libraries(topological_inference ${CGAL_LIBRARY})
	target_link_libraries(morphological_operation_distance_from_point_cloud ${CGAL_LIBRARY})
endif (CGAL_FOUND)

if (TBB_FOUND)
	target_link_libraries(topological_inference ${TBB_LIBRARIES})
	target_link_libraries(topological_inference_periodic_domain ${TBB_LIBRARIES})
	target_link_libraries(morphological_operation_on_cubical_complex ${TBB_LIBRARIES})
	target_link_libraries(morphological_operation_distance_from_point_cloud ${TBB_LIBRARIES})
endif(TBB_FOUND)

