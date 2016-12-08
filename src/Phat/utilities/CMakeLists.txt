cmake_minimum_required(VERSION 2.6)
project(Phat_interface_examples)

if(PHAT_FOUND)
  add_executable ( rips_persistence_with_PHAT rips_persistence_with_PHAT.cpp )
  target_link_libraries(rips_persistence_with_PHAT ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  add_executable ( cubical_complex_persistence_with_PHAT cubical_complex_persistence_with_PHAT.cpp )
  target_link_libraries(cubical_complex_persistence_with_PHAT ${Boost_SYSTEM_LIBRARY})

  add_executable ( cubical_complex_with_periodic_boundary_conditions_persistence_with_PHAT cubical_complex_with_periodic_boundary_conditions_persistence_with_PHAT.cpp )
  target_link_libraries(cubical_complex_with_periodic_boundary_conditions_persistence_with_PHAT ${Boost_SYSTEM_LIBRARY})

  #Witness complex & persistence is currently under review. This example will be added as soon as we have the witness complex merged to a trunk.
  #cmake list for witness complex get a bit more complicated. That is why, at the moment, I am dropping it. 
  #add_executable ( Witness_complex_with_PHAT Witness_complex_with_PHAT.cpp )
  #target_link_libraries(Witness_complex_with_PHAT ${Boost_SYSTEM_LIBRARY})

  #For examples that uses alpha complex we need CGAL and EIGEN
  if(CGAL_FOUND)
    if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
      message(STATUS "CGAL version: ${CGAL_VERSION}.")
      find_package(Eigen3 3.1.0)
      if (EIGEN3_FOUND)
        message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
        include( ${EIGEN3_USE_FILE} )
        add_executable ( alpha_complex_3d_persistence_with_PHAT alpha_complex_3d_persistence_with_PHAT.cpp )
        target_link_libraries(alpha_complex_3d_persistence_with_PHAT ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
        add_executable ( alpha_complex_persistence_with_PHAT alpha_complex_persistence_with_PHAT.cpp )  
        target_link_libraries(alpha_complex_persistence_with_PHAT ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY}) 
      else()
        message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
      endif(EIGEN3_FOUND)
    else()
      message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
    endif(NOT CGAL_VERSION VERSION_LESS 4.7.0)
  endif(CGAL_FOUND)

  if (TBB_FOUND)
    target_link_libraries(rips_persistence_with_PHAT ${TBB_LIBRARIES})
    target_link_libraries(cubical_complex_persistence_with_PHAT ${TBB_LIBRARIES})
    target_link_libraries(cubical_complex_with_periodic_boundary_conditions_persistence_with_PHAT ${TBB_LIBRARIES})
    #target_link_libraries(Witness_complex_with_PHAT ${TBB_LIBRARIES})
    target_link_libraries(alpha_complex_3d_persistence_with_PHAT ${TBB_LIBRARIES})
    target_link_libraries(alpha_complex_persistence_with_PHAT ${TBB_LIBRARIES})
  endif(TBB_FOUND)

else()
    message(WARNING "PHAT was not found in your system. The programs that require PHAT will not be compiled.")
endif(PHAT_FOUND)
