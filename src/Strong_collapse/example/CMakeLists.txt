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
  
  add_executable (rips_persistence_with_sc rips_persistence_with_sc.cpp )
  target_link_libraries(rips_persistence_with_sc ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(rips_persistence_with_sc ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(rips_persistence_with_sc ${TBB_LIBRARIES})
  endif()

  add_executable (generate_sphere_product generate_sphere_product.cpp )
  target_link_libraries(generate_sphere_product ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(generate_sphere_product ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(generate_sphere_product ${TBB_LIBRARIES})
  endif()

  add_executable (manualTowerExample manualTowerExample.cpp )
  target_link_libraries(manualTowerExample ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(manualTowerExample ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(manualTowerExample ${TBB_LIBRARIES})
  endif()
  
  add_executable (compute_bottelneck compute_bottelneck.cpp )
  target_link_libraries(compute_bottelneck ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(compute_bottelneck ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(compute_bottelneck ${TBB_LIBRARIES})
  endif()
  
  
	
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)

