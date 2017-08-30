cmake_minimum_required(VERSION 2.6)
project(GUDHIdev)

include(CMakeGUDHIVersion.txt)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")

enable_testing()

# This variable is used by Cython CMakeLists.txt to know its path
set(GUDHI_CYTHON_PATH "src/cython")
# For third parties libraries management - To be done last as CGAL updates CMAKE_MODULE_PATH
include(GUDHI_third_party_libraries)

if(MSVC)
  # Turn off some VC++ warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("++ Debug compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("++ Release compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (DEBUG_TRACES)
  # For programs to be more verbose
  message(STATUS "DEBUG_TRACES are activated")
  add_definitions(-DDEBUG_TRACES)
endif()

# Add your new module in the list, order is not important
include(GUDHI_modules)

add_gudhi_module(common)
add_gudhi_module(Alpha_complex)
add_gudhi_module(Bitmap_cubical_complex)
add_gudhi_module(Bottleneck_distance)
add_gudhi_module(Contraction)
add_gudhi_module(Hasse_complex)
add_gudhi_module(Persistent_cohomology)
add_gudhi_module(Rips_complex)
add_gudhi_module(Simplex_tree)
add_gudhi_module(Skeleton_blocker)
add_gudhi_module(Spatial_searching)
add_gudhi_module(Subsampling)
add_gudhi_module(Tangential_complex)
add_gudhi_module(Witness_complex)

message("++ GUDHI_MODULES list is:\"${GUDHI_MODULES}\"")

# Include module CMake subdirectories
# GUDHI_SUB_DIRECTORIES is managed in CMAKE_MODULE_PATH/GUDHI_modules.cmake
foreach(GUDHI_MODULE ${GUDHI_MODULES})
  foreach(GUDHI_SUB_DIRECTORY ${GUDHI_SUB_DIRECTORIES})
    if(EXISTS ${CMAKE_SOURCE_DIR}/src/${GUDHI_MODULE}/${GUDHI_SUB_DIRECTORY}/CMakeLists.txt)
      add_subdirectory(src/${GUDHI_MODULE}/${GUDHI_SUB_DIRECTORY}/)
    endif()
  endforeach()
endforeach()

add_subdirectory(src/GudhUI)

if (WITH_GUDHI_PYTHON)
  # specific for cython module
  add_subdirectory(${GUDHI_CYTHON_PATH})
endif()

  # Find TBB package for parallel sort - not mandatory, just optional.
  set(TBB_FIND_QUIETLY ON)
  find_package(TBB)
  if (TBB_FOUND)
    include(${TBB_USE_FILE})
    message("TBB found in ${TBB_LIBRARY_DIRS}")
    add_definitions(-DGUDHI_USE_TBB)
  endif()

  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
    include( ${EIGEN3_USE_FILE} )
    #include_directories (BEFORE "../../include")
  endif (EIGEN3_FOUND)

  # Required programs for unitary tests purpose
  FIND_PROGRAM( GCOVR_PATH gcovr )
  if (GCOVR_PATH)
    message("gcovr found in ${GCOVR_PATH}")
  endif()
  # Required programs for unitary tests purpose
  FIND_PROGRAM( GPROF_PATH gprof )
  if (GPROF_PATH)
    message("gprof found in ${GPROF_PATH}")
  endif()
  FIND_PROGRAM( DIFF_PATH diff )
  if (DIFF_PATH)
    message("diff found in ${DIFF_PATH}")
  endif()

  # BOOST ISSUE result_of vs C++11
  add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
  # BOOST ISSUE with Libraries name resolution under Windows
  add_definitions(-DBOOST_ALL_NO_LIB)
  # problem with Visual Studio link on Boost program_options
  add_definitions( -DBOOST_ALL_DYN_LINK )

  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})                      

  message(STATUS "boost include dirs:" ${Boost_INCLUDE_DIRS})
  message(STATUS "boost library dirs:" ${Boost_LIBRARY_DIRS})

  if (DEBUG_TRACES)
    # For programs to be more verbose
    message(STATUS "DEBUG_TRACES are activated")
    add_definitions(-DDEBUG_TRACES)
  endif()

  include_directories(src/common/include/)
  include_directories(src/Alpha_complex/include/)
  include_directories(src/Bitmap_cubical_complex/include/)
  include_directories(src/Bottleneck_distance/include/)
  include_directories(src/Contraction/include/)
  include_directories(src/Hasse_complex/include/)
  include_directories(src/Persistent_cohomology/include/)
  include_directories(src/Rips_complex/include/)
  include_directories(src/Simplex_tree/include/)
  include_directories(src/Simplex_array_list/include/)
  include_directories(src/Skeleton_blocker/include/)
  include_directories(src/Spatial_searching/include/)
  include_directories(src/Subsampling/include/)
  include_directories(src/Tangential_complex/include/)
  include_directories(src/Witness_complex/include/)

  add_subdirectory(src/common/example)
  add_subdirectory(src/common/test)
  add_subdirectory(src/Simplex_tree/test)
  add_subdirectory(src/Simplex_tree/example)
  add_subdirectory(src/Simplex_array_list/test)
  add_subdirectory(src/Simplex_array_list/example)
  add_subdirectory(src/Persistent_cohomology/test)
  add_subdirectory(src/Persistent_cohomology/example)
  add_subdirectory(src/Persistent_cohomology/benchmark)
  add_subdirectory(src/Skeleton_blocker/test) 
  add_subdirectory(src/Skeleton_blocker/example) 
  add_subdirectory(src/Contraction/example)
  add_subdirectory(src/Witness_complex/test)
  add_subdirectory(src/Witness_complex/example)
  add_subdirectory(src/Bitmap_cubical_complex/test)
  add_subdirectory(src/Bitmap_cubical_complex/example)
  add_subdirectory(src/Alpha_complex/example)
  add_subdirectory(src/Alpha_complex/test)
  add_subdirectory(src/Spatial_searching/example)
  add_subdirectory(src/Spatial_searching/test)
  add_subdirectory(src/Subsampling/example)
  add_subdirectory(src/Subsampling/test)
  add_subdirectory(src/Tangential_complex/example)
  add_subdirectory(src/Tangential_complex/test)
  add_subdirectory(src/Tangential_complex/benchmark)
  add_subdirectory(src/Bottleneck_distance/example)
  add_subdirectory(src/Bottleneck_distance/test)
  add_subdirectory(src/Bottleneck_distance/benchmark)
  add_subdirectory(src/Rips_complex/example)
  add_subdirectory(src/Rips_complex/test)

  # data points generator
  add_subdirectory(data/points/generator)

  # Please let GudhUI in last compilation position as QT is known to modify CMAKE_CXX_FLAGS
  # GudhUI
  add_subdirectory(src/GudhUI)

endif()  


