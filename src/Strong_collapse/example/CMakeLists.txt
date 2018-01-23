cmake_minimum_required(VERSION 2.6)
project(Strong_collapse)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  
  add_executable (sparseCollapseRips sparseCollapseRips.cpp )
  target_link_libraries(sparseCollapseRips ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(sparseCollapseRips ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(sparseCollapseRips ${TBB_LIBRARIES})
  endif()
  
  add_executable (manualSimplicialComplexColl manualSimplicialComplexColl.cpp )
  target_link_libraries(manualSimplicialComplexColl ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(manualSimplicialComplexColl ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(manualSimplicialComplexColl ${TBB_LIBRARIES})
  endif()

  add_executable (alphaSparseCollapse alphaSparseCollapse.cpp )
  target_link_libraries(alphaSparseCollapse ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(alphaSparseCollapse ${Boost_PROGRAM_OPTIONS_LIBRARY})		
  if (TBB_FOUND)
   target_link_libraries(alphaSparseCollapse ${TBB_LIBRARIES})
  endif()	
  
  add_executable (rips_for_sophia rips_for_sophia.cpp )
  target_link_libraries(rips_for_sophia ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(rips_for_sophia ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(rips_for_sophia ${TBB_LIBRARIES})
  endif()

  add_executable (generate_sphere_product generate_sphere_product.cpp )
  target_link_libraries(generate_sphere_product ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(generate_sphere_product ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(generate_sphere_product ${TBB_LIBRARIES})
  endif()		
 		
	
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)

