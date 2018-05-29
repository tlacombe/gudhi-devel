project(Strong_collapse_examples)

add_executable ( Strong_collapse_example_mini_simplex_tree mini_simplex_tree.cpp )
if (TBB_FOUND)
  target_link_libraries(Strong_collapse_example_mini_simplex_tree ${TBB_LIBRARIES})
endif()
