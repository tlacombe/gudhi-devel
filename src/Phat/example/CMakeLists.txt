cmake_minimum_required(VERSION 2.6)
project(GUDHIBitmap)

add_executable ( rips_persistence_with_phat rips_persistence_with_phat.cpp )
target_link_libraries(rips_persistence_with_phat ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
#add_test(Random_bitmap_cubical_complex ${CMAKE_CURRENT_BINARY_DIR}/Random_bitmap_cubical_complex 2 100 100)	

add_executable ( cubical_complex_persistence_with_phat cubical_complex_persistence_with_phat.cpp )
target_link_libraries(cubical_complex_persistence_with_phat ${Boost_SYSTEM_LIBRARY})
#add_test(Random_bitmap_cubical_complex ${CMAKE_CURRENT_BINARY_DIR}/Random_bitmap_cubical_complex 2 100 100)


