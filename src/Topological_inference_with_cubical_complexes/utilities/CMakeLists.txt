cmake_minimum_required(VERSION 2.6)
project(GUDHITopologicalInferenceWithCubicalComplexes)


#if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
   add_executable ( topological_inference topological_inference.cpp )
   target_link_libraries(topological_inference ${Boost_SYSTEM_LIBRARY})
   
   #add_executable ( distance_from_collection_of_cubes distance_from_collection_of_cubes.cpp )	
   #target_link_libraries(distance_from_collection_of_cubes ${Boost_SYSTEM_LIBRARY})

if ( CGAL_WITH_EIGEN3_VERSION VERSION_MORE 4.6.0 )
	target_link_libraries(topological_inference ${TBB_LIBRARIES})	
endif( CGAL_WITH_EIGEN3_VERSION VERSION_MORE 4.6.0 )

   if (TBB_FOUND)
      target_link_libraries(topological_inference ${TBB_LIBRARIES})
      #target_link_libraries(distance_from_collection_of_cubes ${TBB_LIBRARIES})
   endif()
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
