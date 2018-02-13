cmake_minimum_required(VERSION 2.6)
project(Bitmap_cubical_complex_examples)

add_executable ( Random_bitmap_cubical_complex Random_bitmap_cubical_complex.cpp )
if (TBB_FOUND)
  target_link_libraries(Random_bitmap_cubical_complex ${TBB_LIBRARIES})
endif()
add_test(NAME Bitmap_cubical_complex_example_random COMMAND $<TARGET_FILE:Random_bitmap_cubical_complex>
    "2" "100" "100")


add_executable ( setting_up_filtration_of_individual_cells_example setting_up_filtration_of_individual_cells_example.cpp )
if (TBB_FOUND)
  target_link_libraries(setting_up_filtration_of_individual_cells_example ${TBB_LIBRARIES})
endif()
add_test(NAME setting_up_filtration_of_individual_cells_example COMMAND $<TARGET_FILE:setting_up_filtration_of_individual_cells_example>)

install(TARGETS Random_bitmap_cubical_complex DESTINATION bin)
