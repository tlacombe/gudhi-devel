cmake_minimum_required(VERSION 2.6)
project(GUDHIExPersCohom)

find_package(Boost 1.45.0 COMPONENTS system unit_test_framework program_options)
#TODO: find_package(gmp REQUIRED COMPONENTS)

if(MSVC)
	# Turn off some VC++ warnings
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -Wall -lgmpxx -lgmp -lboost_unit_test_framework -lboost_system")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

if(NOT Boost_FOUND)
   message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
   message("Boost_include = ${Boost_INCLUDE_DIRS}")
   message("Boost_lib = ${Boost_LIBRARIES}")

   include_directories(${Boost_INCLUDE_DIRS})
   include_directories(BEFORE "../../common" "../../common/include/gudhi" "../include/gudhi" "../../Simplex_tree/include/gudhi" "../../Hasse_complex/include/gudhi" )
   add_executable ( rips_persistence rips_persistence.cpp )
   target_link_libraries(rips_persistence ${Boost_LIBRARIES})

   if(gmp_FOUND)
      add_executable ( rips_multifield_persistence rips_multifield_persistence.cpp )
      target_link_libraries(rips_multifield_persistence ${Boost_LIBRARIES})

      add_executable ( performance_rips_persistence performance_rips_persistence.cpp )
      target_link_libraries(performance_rips_persistence ${Boost_LIBRARIES})
   endif()
endif()  
