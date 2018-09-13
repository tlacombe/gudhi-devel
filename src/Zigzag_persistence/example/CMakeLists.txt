project(Simplex_tree_examples)

add_executable ( tmp_zz_points tmp_zz_points.cpp )
if (TBB_FOUND)
  target_link_libraries(tmp_zz_points ${TBB_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
else()
  target_link_libraries(tmp_zz_points ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()

add_executable ( Zigzag_persistence_simple_filtration zigzag_persistence_simple_filtration.cpp )
if (TBB_FOUND)
  target_link_libraries(Zigzag_persistence_simple_filtration ${TBB_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
else()
  target_link_libraries(Zigzag_persistence_simple_filtration ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()



add_executable ( Zigzag_persistence_flag_zigzag_points zigzag_persistence_flag_zigzag_points.cpp )
if (TBB_FOUND)
  target_link_libraries(Zigzag_persistence_flag_zigzag_points ${TBB_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
else()
  target_link_libraries(Zigzag_persistence_flag_zigzag_points ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()

