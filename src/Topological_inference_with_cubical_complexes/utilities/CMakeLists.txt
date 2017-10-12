cmake_minimum_required(VERSION 2.6)
project(GUDHITopologicalInferenceWithCubicalComplexes)

add_executable ( topological_inference topological_inference.cpp )
target_link_libraries(topological_inference ${Boost_SYSTEM_LIBRARY})

add_executable ( topological_inference_periodic_domain topological_inference_periodic_domain.cpp )
target_link_libraries(topological_inference_periodic_domain ${Boost_SYSTEM_LIBRARY})
   
#add_executable ( distance_from_collection_of_cubes distance_from_collection_of_cubes.cpp )	
#target_link_libraries(distance_from_collection_of_cubes ${Boost_SYSTEM_LIBRARY})

if (CGAL_FOUND)
	target_link_libraries(topological_inference ${CGAL_LIBRARY})	
endif (CGAL_FOUND)

if (TBB_FOUND)
	target_link_libraries(topological_inference ${TBB_LIBRARIES})
	target_link_libraries(topological_inference_periodic_domain ${TBB_LIBRARIES})
	#target_link_libraries(distance_from_collection_of_cubes ${TBB_LIBRARIES})
endif(TBB_FOUND)

