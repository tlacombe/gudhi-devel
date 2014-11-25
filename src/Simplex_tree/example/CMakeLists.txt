cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeFromFile)

set(CMAKE_PREFIX_PATH "../../cmake/modules/")
set(CMAKE_MODULE_PATH "../../cmake/modules/")
message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

find_package(Boost 1.45.0 COMPONENTS system thread)
find_package(CGAL)
find_package(GMPXX)
find_package(GMP)

if(MSVC)
	# Turn off some VC++ warnings
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
    add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -Werror -lboost_system")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

if(NOT Boost_FOUND)
   message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
   message("Boost_include = ${Boost_INCLUDE_DIRS}")
   message("Boost_lib = ${Boost_LIBRARIES}")

   include_directories(${Boost_INCLUDE_DIRS})
   include_directories(BEFORE "../../common/include/" "../include/" )
   add_executable ( simplex_tree_from_file simplex_tree_from_file.cpp )
   target_link_libraries(simplex_tree_from_file ${Boost_LIBRARIES})
   
   # An example with Simplex-tree using CGAL alpha_shapes_3
   if(GMPXX_FOUND AND GMP_FOUND AND CGAL_FOUND)
      message("CGAL_include = ${CGAL_INCLUDE_DIRS}")
      message("CGAL_lib = ${CGAL_LIBRARIES_DIR}")

      include_directories(${Boost_INCLUDE_DIRS})
      add_executable ( alpha_shapes_3_from_simplex_tree alpha_shapes_3_from_simplex_tree.cpp )
      target_link_libraries(alpha_shapes_3_from_simplex_tree ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY})
   endif()  
   
endif()  
